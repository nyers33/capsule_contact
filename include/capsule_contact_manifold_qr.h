#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

template<typename T> using unsized_raw_array = T[];

// manifold data structure & collision feature uses the pattern
// from Box2D-Lite by Eric Catto

// capsule collision feature numbering
//
//         ^ y
//         |
//         1
//      /-----\
//     /       \
//    /         \
//   |           |
//   |           |
//   |     0     | --> x
//   |           |
//   |           |
//    \         /
//     \       /
//      \-----/
//         -1

enum CollisionFeature
{
	NEG_Y = -1,
	BODY,
	POS_Y
};

union FeaturePair
{
	struct Feature
	{
		char inFeature;
		char outFeature;
	} e;
	int value;
};

template <typename Real, int N>
struct Contact
{
	Contact() : separation(static_cast<Real>(0)), feature{ 0 } {}

	Eigen::Matrix<Real, N, 1> position;
	Eigen::Matrix<Real, N, 1> normal;
	Real separation;
	FeaturePair feature;
};

template <typename Real, int N>
struct Capsule
{
	Capsule() : height(static_cast<Real>(1)), rad(static_cast<Real>(0.5)) {}
	Capsule(const Eigen::Matrix<Real, N+1, N+1>& pose, Real height, Real rad) : pose(pose), height(height), rad(rad) {}

	Eigen::Matrix<Real, N+1, N+1> pose;
	Real height, rad;

	auto translation() const
	{
		return pose.col(N).block<N, 1>(0, 0);
	}

	auto rotation() const
	{
		return pose.block<N, N>(0, 0);
	}
};

// QR Decomposition with the Gram-Schmidt Procedure from RPubs
template <typename Real, int N>
void qrDecompGS(const Eigen::Matrix<Real, N, 2>& inMat, Eigen::Matrix<Real, N, 2>& qMat, Eigen::Matrix<Real, 2, 2>& rMat)
{
	auto zeroVec = Eigen::Matrix<Real, N, 1>::Zero();
	const Real eps = static_cast<Real>(1e-4);

	auto a1 = inMat.col(0);
	auto a2 = inMat.col(1);

	auto v1 = a1;
	Real v1Norm = v1.norm();
	auto e1 = (v1Norm > eps) ? v1.normalized() : zeroVec;
	auto v2 = a2 - a2.dot(e1) * e1;
	Real v2Norm = v2.norm();
	auto e2 = (v2Norm > eps) ? v2.normalized() : zeroVec;

	qMat.col(0) = e1;
	qMat.col(1) = e2;

	rMat(0, 0) = a1.dot(e1);
	rMat(1, 0) = static_cast<Real>(0);
	rMat(0, 1) = a2.dot(e1);
	rMat(1, 1) = a2.dot(e2);
}

// Real-Time Collision Detection by Christer Ericson
template <typename Real>
Real sqDistPointSeg(Eigen::Matrix<Real, 2, 1> a, Eigen::Matrix<Real, 2, 1> b, Eigen::Matrix<Real, 2, 1> c)
{
	typedef Eigen::Matrix<Real, 2, 1> Vec2;

	Vec2 ab = b - a;
	Vec2 ac = c - a;
	Vec2 bc = c - b;
	Real e = ac.dot(ab);
	Real f = ab.dot(ab);
	if (e < static_cast<Real>(0))
	{
		return ac.dot(ac);
	}
	else if (e >= f)
	{
		return bc.dot(bc);
	}
	else
	{
		return ac.dot(ac) - e * e / f;
	}
}

// Real-Time Collision Detection by Christer Ericson
template <typename Real>
Eigen::Matrix<Real, 2, 1> closestPtPointSeg(Eigen::Matrix<Real, 2, 1> a, Eigen::Matrix<Real, 2, 1> b, Eigen::Matrix<Real, 2, 1> c)
{
	typedef Eigen::Matrix<Real, 2, 1> Vec2;
	const Real eps = static_cast<Real>(1e-4);

	Vec2 ab = b - a;
	Real t = (c - a).dot(ab) / ab.dot(ab);
	if (ab.dot(ab) < eps)
		t = static_cast<Real>(0);
	t = std::clamp(t, static_cast<Real>(0), static_cast<Real>(1));
	return a + t * ab;
}

// find the point on one of the parallelogram's edges that is closest to given point
// if given point is inside the parallelogram, return given point
template <typename Real>
Eigen::Matrix<Real, 2, 1> parallelogramContainsPt(Eigen::Matrix<Real, 2, 1> o, Eigen::Matrix<Real, 2, 1> a, Eigen::Matrix<Real, 2, 1> b, Eigen::Matrix<Real, 2, 1> pt)
{
	typedef Eigen::Matrix<Real, 2, 1> Vec2;
	const Real eps = static_cast<Real>(1e-4);

	Vec2 p = pt - o;
	Real denomAbs = std::abs(a.x() * b.y() - a.y() * b.x());

	Real mu, lambda;

	mu = (p.x() * b.y() - p.y() * b.x()) / (a.x() * b.y() - a.y() * b.x());
	lambda = (p.x() * a.y() - p.y() * a.x()) / (a.y() * b.x() - a.x() * b.y());

	if (denomAbs > eps && (static_cast<Real>(0) <= mu && mu <= static_cast<Real>(1) && static_cast<Real>(0) <= lambda && lambda <= static_cast<Real>(1)))
	{
		// inside
		return pt;
	}
	else
	{
		// outside
		Vec2 edgeList[4][2] = { { Vec2(0, 0), Vec2(a.x(), a.y()) }, { Vec2(a.x(), a.y()), Vec2(a.x(), a.y()) + Vec2(b.x(), b.y()) }, { Vec2(a.x(), a.y()) + Vec2(b.x(), b.y()), Vec2(b.x(), b.y()) }, { Vec2(b.x(), b.y()), Vec2(0, 0) } };
		Real minDistToEdge = std::numeric_limits<Real>::max();
		int iEdge = -1;
		for (int i = 0; i < 4; i++)
		{
			Real distToEdge = sqDistPointSeg(edgeList[i][0], edgeList[i][1], p);
			if (distToEdge <= minDistToEdge)
			{
				minDistToEdge = distToEdge;
				iEdge = i;
			}
		}

		return closestPtPointSeg(edgeList[iEdge][0], edgeList[iEdge][1], p) + o;
	}
}

// Efficient Calculation of Minimum Distance Between Capsules and Its Use in Robotics
template <typename Real, int N>
Real distanceQR(const Eigen::Matrix<Real, N, 1> bodyA[2], const Eigen::Matrix<Real, N, 1> bodyB[2], Eigen::Matrix<Real, N, 2>* qMat_out = nullptr, Eigen::Matrix<Real, 2, 2>* rMat_out = nullptr, Eigen::Matrix<Real, 2, 1>* opt_out = nullptr)
{
	typedef Eigen::Matrix<Real, 2, 1> Vec2;

	const auto& p1 = bodyA[0]; // seg0p0
	const auto& u1 = bodyA[1]; // seg0p1

	const auto& p2 = bodyB[0]; // seg1p0
	const auto& u2 = bodyB[1]; // seg1p1

	auto s1 = u1 - p1;
	auto s2 = u2 - p2;

	Eigen::Matrix<Real, N, 2> inMat, qMat;
	Eigen::Matrix<Real, 2, 2> rMat;
	inMat.col(0) = s2;
	inMat.col(1) = -s1;

	qrDecompGS(inMat, qMat, rMat);
	if (qMat_out)
	{
		*qMat_out = qMat;
	}
	if (rMat_out)
	{
		*rMat_out = rMat;
	}

	Vec2 parallelogramOrigin = rMat * Vec2(0, 0) + qMat.transpose() * (p2 - p1);
	Vec2 parallelogramSideA = rMat * Vec2(1, 0) + qMat.transpose() * (p2 - p1) - parallelogramOrigin;
	Vec2 parallelogramSideB = rMat * Vec2(0, 1) + qMat.transpose() * (p2 - p1) - parallelogramOrigin;

	Vec2 opt = parallelogramContainsPt(parallelogramOrigin, parallelogramSideA, parallelogramSideB, Vec2(0, 0));
	if (opt_out)
	{
		*opt_out = opt;
	}

	// sqDist
	return opt.dot(opt) + (p2 - p1).dot(p2 - p1) - (p2 - p1).dot(qMat * qMat.transpose() * (p2 - p1));
}

template <typename Real, int N>
int manifoldQR(Contact<Real, N>* contacts, const Capsule<Real, N>& bodyA, const Capsule<Real, N>& bodyB)
{
	typedef Eigen::Matrix<Real, 2, 1> Vec2;
	const auto zeroVec = Eigen::Matrix<Real, N, 1>::Zero();
	const Real eps = static_cast<Real>(1e-4);
	
	auto p1 = bodyA.translation() - static_cast<Real>(0.5) * bodyA.height * bodyA.rotation().col(1); // seg0p0
	auto u1 = bodyA.translation() + static_cast<Real>(0.5) * bodyA.height * bodyA.rotation().col(1); // seg0p1

	auto p2 = bodyB.translation() - static_cast<Real>(0.5) * bodyB.height * bodyB.rotation().col(1); // seg1p0
	auto u2 = bodyB.translation() + static_cast<Real>(0.5) * bodyB.height * bodyB.rotation().col(1); // seg1p1

	auto s1 = u1 - p1;
	auto s2 = u2 - p2;

	Real r1 = bodyA.rad;
	Real r2 = bodyB.rad;

	auto&& p1u1 = unsized_raw_array< Eigen::Matrix<Real, N, 1> >{ p1, u1 };
	auto&& p2u2 = unsized_raw_array< Eigen::Matrix<Real, N, 1> >{ p2, u2 };
	
	// 3D - Mat32;
	// 2D - Mat22;
	auto qMat = Eigen::Matrix<Real, N, 2>();
	Eigen::Matrix<Real, 2, 2> rMat;
	Vec2 opt;
	
	Real sqDist = distanceQR(p1u1, p2u2, &qMat, &rMat, &opt);

	if ((bodyA.rad + bodyB.rad) * (bodyA.rad + bodyB.rad) < sqDist)
	{
		return 0;
	}

	Vec2 rx = opt - qMat.transpose() * (p2 - p1);

	Real parallelTest = std::abs(std::abs(s1.dot(s2)) - s1.norm() * s2.norm());

	// line vs line nearest in point - not parallel
	if (std::abs(rMat.col(0).x()) > eps && std::abs(rMat.col(1).y()) > eps && parallelTest > eps)
	{
		Real x1 = rx.y() / rMat.col(1).y();
		Real x2 = (rMat.col(1).y() * rx.x() - rMat.col(1).x() * rx.y()) / (rMat.col(0).x() * rMat.col(1).y());

		if (x2 < static_cast<Real>(0))
			x2 = static_cast<Real>(0);

		FeaturePair fp;
		if (x1 <= static_cast<Real>(0) + eps)
		{
			fp.e.inFeature = NEG_Y;
		}
		else if (x1 >= static_cast<Real>(1) - eps)
		{
			fp.e.inFeature = POS_Y;
		}
		else
		{
			fp.e.inFeature = BODY;
		}
		if (x2 <= static_cast<Real>(0) + eps)
		{
			fp.e.outFeature = NEG_Y;
		}
		else if (x2 >= static_cast<Real>(1) - eps)
		{
			fp.e.outFeature = POS_Y;
		}
		else
		{
			fp.e.outFeature = BODY;
		}

		auto nContact = ((p2 + x2 * s2) - (p1 + x1 * s1));
		Real nContactNorm = nContact.norm();
		if (nContactNorm > eps)
		{
			contacts[0].normal = nContact.normalized();
		}
		else
		{
			contacts[0].normal = zeroVec;
		}

		contacts[0].separation = sqrtf(sqDist) - r1 - r2;
		contacts[0].position = static_cast<Real>(0.5) * (((p1 + x1 * s1) + r1 * contacts[0].normal) + ((p2 + x2 * s2) - r2 * contacts[0].normal));
		contacts[0].feature = fp;

		return 1;
	}
	
	Real s1LenSqd = s1.dot(s1);
	Real s2LenSqd = s2.dot(s2);
	
	// point vs point
	if (s1LenSqd < eps && s2LenSqd < eps)
	{
		FeaturePair fp;
		fp.e.inFeature = BODY;
		fp.e.outFeature = BODY;

		auto nContact = (p2 - p1);
		Real nContactNorm = nContact.norm();
		if (nContactNorm > eps)
		{
			contacts[0].normal = nContact.normalized();
		}
		else
		{
			contacts[0].normal = zeroVec;
		}

		contacts[0].separation = sqrtf(sqDist) - r1 - r2;
		contacts[0].position = static_cast<Real>(0.5) * ((p1 + r1 * contacts[0].normal) + (p2 - r2 * contacts[0].normal));
		contacts[0].feature = fp;

		return 1;
	}
	
	// line vs line - parallel lines
	if (s1LenSqd * s2LenSqd > eps)
	{
		Real x1[2];
		Real x2[2];
		FeaturePair fp;

		if (s1LenSqd > s2LenSqd)
		{
			x1[0] = std::min(std::max(((static_cast<Real>(1) / sqrtf(s1LenSqd)) * s1).dot(p2 - p1) / sqrtf(s1LenSqd), static_cast<Real>(0)), static_cast<Real>(1));
			x1[1] = std::min(std::max(((static_cast<Real>(1) / sqrtf(s1LenSqd)) * s1).dot(u2 - p1) / sqrtf(s1LenSqd), static_cast<Real>(0)), static_cast<Real>(1));

			x2[0] = std::min(std::max(((static_cast<Real>(1) / sqrtf(s2LenSqd)) * s2).dot(p1 - p2 + x1[0] * s1) / sqrtf(s2LenSqd), static_cast<Real>(0)), static_cast<Real>(1));
			x2[1] = std::min(std::max(((static_cast<Real>(1) / sqrtf(s2LenSqd)) * s2).dot(p1 - p2 + x1[1] * s1) / sqrtf(s2LenSqd), static_cast<Real>(0)), static_cast<Real>(1));

		}
		else
		{
			x2[0] = std::min(std::max(((static_cast<Real>(1) / sqrtf(s2LenSqd)) * s2).dot(p1 - p2) / sqrtf(s2LenSqd), static_cast<Real>(0)), static_cast<Real>(1));
			x2[1] = std::min(std::max(((static_cast<Real>(1) / sqrtf(s2LenSqd)) * s2).dot(u1 - p2) / sqrtf(s2LenSqd), static_cast<Real>(0)), static_cast<Real>(1));

			x1[0] = std::min(std::max(((static_cast<Real>(1) / sqrtf(s1LenSqd)) * s1).dot(p2 - p1 + x2[0] * s2) / sqrtf(s1LenSqd), static_cast<Real>(0)), static_cast<Real>(1));
			x1[1] = std::min(std::max(((static_cast<Real>(1) / sqrtf(s1LenSqd)) * s1).dot(p2 - p1 + x2[1] * s2) / sqrtf(s1LenSqd), static_cast<Real>(0)), static_cast<Real>(1));
		}

		auto nContact = (p2 + static_cast<Real>(0.5) * (x2[0] + x2[1]) * s2) - (p1 + static_cast<Real>(0.5) * (x1[0] + x1[1]) * s1);
		Real nContactNorm = nContact.norm();
		if (nContactNorm > eps)
		{
			contacts[0].normal = nContact.normalized();
		}
		else
		{
			contacts[0].normal = zeroVec;
		}

		for (int i = 0; i < 2; ++i)
		{
			contacts[i].separation = sqrtf(sqDist) - r1 - r2;
			contacts[i].position = (p1 + x1[i] * s1) + static_cast<Real>(0.5) * (r1 - r2 + sqrtf(sqDist)) * contacts[i].normal;
			fp.e.inFeature = BODY;
			fp.e.outFeature = BODY;
			contacts[i].feature = fp;
		}

		// capsule cap check
		if (
			((x1[0] <= static_cast<Real>(0) + eps && x1[1] <= static_cast<Real>(0) + eps) || (x1[0] >= static_cast<Real>(1) - eps && x1[1] >= static_cast<Real>(1) - eps))
			&&
			((x2[0] <= static_cast<Real>(0) + eps && x2[1] <= static_cast<Real>(0) + eps) || (x2[0] >= static_cast<Real>(1) - eps && x2[1] >= static_cast<Real>(1) - eps))
			)
		{
			if (x1[0] <= static_cast<Real>(0) + eps)
			{
				fp.e.inFeature = NEG_Y;
			}
			if (x1[0] >= static_cast<Real>(1) - eps)
			{
				fp.e.inFeature = POS_Y;
			}
			if (x2[0] <= static_cast<Real>(0) + eps)
			{
				fp.e.outFeature = NEG_Y;
			}
			if (x2[0] >= static_cast<Real>(1) - eps)
			{
				fp.e.outFeature = POS_Y;
			}
			contacts[0].feature = fp;

			return 1;
		}
		else
		{
			return 2;
		}
	}
	
	// line vs point
	{
		Real x1[2];
		Real x2[2];
		FeaturePair fp;

		if (s1LenSqd > s2LenSqd)
		{
			x1[1] = std::min(std::max(((static_cast<Real>(1) / sqrtf(s1LenSqd)) * s1).dot(u2 - p1) / sqrtf(s1LenSqd), static_cast<Real>(0)), static_cast<Real>(1));
			x1[0] = std::min(std::max(((static_cast<Real>(1) / sqrtf(s1LenSqd)) * s1).dot(p2 - p1) / sqrtf(s1LenSqd), static_cast<Real>(0)), static_cast<Real>(1));

			x2[0] = static_cast<Real>(0);
			x2[1] = static_cast<Real>(1);

			if (x1[0] == static_cast<Real>(0))
			{
				fp.e.inFeature = NEG_Y;
			}
			else if (x1[0] == static_cast<Real>(1))
			{
				fp.e.inFeature = POS_Y;
			}
			else
			{
				fp.e.inFeature = BODY;
			}
			fp.e.outFeature = BODY;
		}
		else
		{
			x1[0] = static_cast<Real>(0);
			x1[1] = static_cast<Real>(1);

			x2[0] = std::min(std::max(((static_cast<Real>(1) / sqrtf(s2LenSqd)) * s2).dot(p1 - p2) / sqrtf(s2LenSqd), static_cast<Real>(0)), static_cast<Real>(1));
			x2[1] = std::min(std::max(((static_cast<Real>(1) / sqrtf(s2LenSqd)) * s2).dot(u1 - p2) / sqrtf(s2LenSqd), static_cast<Real>(0)), static_cast<Real>(1));

			if (x2[0] == static_cast<Real>(0))
			{
				fp.e.outFeature = NEG_Y;
			}
			else if (x2[0] == static_cast<Real>(1))
			{
				fp.e.outFeature = POS_Y;
			}
			else
			{
				fp.e.outFeature = BODY;
			}
			fp.e.inFeature = BODY;
		}
		
		auto nContact = ((p2 + x2[0] * s2) - (p1 + x1[0] * s1));
		Real nContactNorm = nContact.norm();
		if (nContactNorm > eps)
		{
			contacts[0].normal = nContact.normalized();
		}
		else
		{
			contacts[0].normal = zeroVec;
		}
		
		contacts[0].separation = sqrtf(sqDist) - r1 - r2;
		contacts[0].position = (p1 + x1[0] * s1) + static_cast<Real>(0.5) * (r1 - r2 + sqrtf(sqDist)) * contacts[0].normal;
		contacts[0].feature = fp;
		
		return 1;
	}

	return 0;
}

// 2D rotation matrix around the origin
template <typename Real>
Eigen::Matrix<Real, 2, 2> rotation2D(Real angle)
{
	Eigen::Matrix<Real, 2, 2> rotMat = Eigen::Matrix<Real, 2, 2>();
	rotMat(0, 0) = std::cos(angle);
	rotMat(0, 1) = -std::sin(angle);
	rotMat(1, 0) = std::sin(angle);
	rotMat(1, 1) = std::cos(angle);

	return rotMat;
}

// Euler rotation matrix
// rotation order XYZ
template <typename Real>
Eigen::Matrix<Real, 3, 3> rotation3D(Real a, Real b, Real c)
{
	Eigen::Matrix<Real, 3, 3> eulerMat = Eigen::Matrix<Real, 3, 3>();
	Real sin_a = std::sin(a);
	Real cos_a = std::cos(a);
	Real sin_b = std::sin(b);
	Real cos_b = std::cos(b);
	Real sin_c = std::sin(c);
	Real cos_c = std::cos(c);
	
	eulerMat(0, 0) = cos_b * cos_c;
	eulerMat(0, 1) = -cos_b * sin_c;
	eulerMat(0, 2) = sin_b;

	eulerMat(1, 0) = cos_c * sin_a * sin_b + cos_a * sin_c;
	eulerMat(1, 1) = cos_a * cos_c - sin_a * sin_b * sin_c;
	eulerMat(1, 2) = -cos_b * sin_a;

	eulerMat(2, 0) = -cos_a * cos_c * sin_b + sin_a * sin_c;
	eulerMat(2, 1) = cos_c * sin_a + cos_a * sin_b * sin_c;
	eulerMat(2, 2) = cos_a * cos_b;

	return eulerMat;
}

// Inigo Quilez
template <typename Real>
Eigen::Matrix<Real, 3, 3> rotationAlign3D(Eigen::Matrix<Real, 3, 1> a, Eigen::Matrix<Real, 3, 1> b)
{
	a.normalize();
	b.normalize();

	const Eigen::Matrix<Real, 3, 1> v = a.cross(b);
	Real c = a.dot(b);
	Real k = static_cast<Real>(1) / (static_cast<Real>(1) + c);

	Eigen::Matrix<Real, 3, 3> eulerMat = Eigen::Matrix<Real, 3, 3>();

	eulerMat <<
		v.x() * v.x() * k + c, v.y() * v.x() * k - v.z(), v.z() * v.x() * k + v.y(),
		v.x() * v.y() * k + v.z(), v.y() * v.y() * k + c, v.z() * v.y() * k - v.x(),
		v.x() * v.z() * k - v.y(), v.y() * v.z() * k + v.x(), v.z() * v.z() * k + c;

	return eulerMat;
}

// 2D homogeneous matrix
template <typename Real>
Eigen::Matrix<Real, 3, 3> homogeneous2D(Eigen::Matrix<Real, 2, 2> linear, Eigen::Matrix<Real, 2, 1> displacement)
{
	Eigen::Matrix<Real, 3, 3> hMat = Eigen::Matrix<Real, 3, 3>::Identity();
	hMat.block<2, 2>(0, 0) = linear;
	hMat.block<2, 1>(0, 2) = displacement;

	return hMat;
}

// 3D homogeneous matrix
template <typename Real>
Eigen::Matrix<Real, 4, 4> homogeneous3D(Eigen::Matrix<Real, 3, 3> linear, Eigen::Matrix<Real, 3, 1> displacement)
{
	Eigen::Matrix<Real, 4, 4> hMat = Eigen::Matrix<Real, 4, 4>::Identity();
	hMat.block<3, 3>(0, 0) = linear;
	hMat.block<3, 1>(0, 3) = displacement;

	return hMat;
}
