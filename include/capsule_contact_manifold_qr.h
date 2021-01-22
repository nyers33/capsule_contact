#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

typedef Eigen::Matrix<float, 2, 2> Mat22;
typedef Eigen::Matrix<float, 3, 3> Mat33;
typedef Eigen::Matrix<float, 2, 1> Vec2;

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

struct Contact
{
	Contact() : separation(0.0f), feature{ 0 } {}

	Vec2 position;
	Vec2 normal;
	float separation;
	FeaturePair feature;
};

struct Capsule
{
	Capsule() : height(1.0f), rad(0.5f) {}
	Capsule(Mat33 pose, float height, float rad) : pose(pose), height(height), rad(rad) {}

	Mat33 pose;
	float height, rad;

	Vec2 translation() const
	{
		return pose.block<2, 1>(0, 2);
	}

	Mat22 rotation() const
	{
		return pose.block<2, 2>(0, 0);
	}
};

// QR Decomposition with the Gram-Schmidt Procedure from RPubs
void qrDecompGS(const Mat22& inMat, Mat22& qMat, Mat22& rMat)
{
	Vec2 a1 = inMat.col(0);
	Vec2 a2 = inMat.col(1);

	Vec2 v1 = a1;
	float v1Norm = v1.norm();
	Vec2 e1 = (v1Norm > 1e-4) ? v1.normalized() : Vec2(0.0f, 0.0f);
	Vec2 v2 = a2 - a2.dot(e1) * e1;
	float v2Norm = v2.norm();
	Vec2 e2 = (v2Norm > 1e-4) ? v2.normalized() : Vec2(0.0f, 0.0f);

	qMat.col(0) = e1;
	qMat.col(1) = e2;

	rMat(0, 0) = a1.dot(e1);
	rMat(1, 0) = 0.0f;
	rMat(0, 1) = a2.dot(e1);
	rMat(1, 1) = a2.dot(e2);
}

// Real-Time Collision Detection by Christer Ericson
float sqDistPointSeg(Vec2 a, Vec2 b, Vec2 c)
{
	Vec2 ab = b - a;
	Vec2 ac = c - a;
	Vec2 bc = c - b;
	float e = ac.dot(ab);
	float f = ab.dot(ab);
	if (e < 0.0f)
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
Vec2 closestPtPointSeg(Vec2 a, Vec2 b, Vec2 c)
{
	Vec2 ab = b - a;
	float t = (c - a).dot(ab) / ab.dot(ab);
	if (ab.dot(ab) < 1e-4)
		t = 0.0f;
	t = std::clamp(t, 0.0f, 1.0f);
	return a + t * ab;
}

// find the point on one of the parallelogram's edges that is closest to given point
// if given point is inside the parallelogram, return given point
Vec2 parallelogramContainsPt(Vec2 o, Vec2 a, Vec2 b, Vec2 pt)
{
	Vec2 p = pt - o;
	float denomAbs = std::abs(a.x() * b.y() - a.y() * b.x());

	float mu, lambda;

	mu = (p.x() * b.y() - p.y() * b.x()) / (a.x() * b.y() - a.y() * b.x());
	lambda = (p.x() * a.y() - p.y() * a.x()) / (a.y() * b.x() - a.x() * b.y());

	if (denomAbs > 1e-4 && (0.0f <= mu && mu <= 1.0f && 0.0f <= lambda && lambda <= 1.0f))
	{
		// inside
		return pt;
	}
	else
	{
		// outside
		Vec2 edgeList[4][2] = { { Vec2(0.0f, 0.0f), Vec2(a.x(), a.y()) }, { Vec2(a.x(), a.y()), Vec2(a.x(), a.y()) + Vec2(b.x(), b.y()) }, { Vec2(a.x(), a.y()) + Vec2(b.x(), b.y()), Vec2(b.x(), b.y()) }, { Vec2(b.x(), b.y()), Vec2(0.0f, 0.0f) } };
		float minDistToEdge = std::numeric_limits<float>::max();
		int iEdge = -1;
		for (int i = 0; i < 4; i++)
		{
			float distToEdge = sqDistPointSeg(edgeList[i][0], edgeList[i][1], p);
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
float distanceQR(const Vec2 bodyA[2], const Vec2 bodyB[2], Mat22* qMat_out = nullptr, Mat22* rMat_out = nullptr, Vec2* opt_out = nullptr)
{
	const Vec2& p1 = bodyA[0]; // seg0p0
	const Vec2& u1 = bodyA[1]; // seg0p1

	const Vec2& p2 = bodyB[0]; // seg1p0
	const Vec2& u2 = bodyB[1]; // seg1p1

	Vec2 s1 = u1 - p1;
	Vec2 s2 = u2 - p2;

	Mat22 inMat, qMat, rMat;
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

	Vec2 parallelogramOrigin = rMat * Vec2(0.0f, 0.0f) + qMat.transpose() * (p2 - p1);
	Vec2 parallelogramSideA = rMat * Vec2(1.0f, 0.0f) + qMat.transpose() * (p2 - p1) - parallelogramOrigin;
	Vec2 parallelogramSideB = rMat * Vec2(0.0f, 1.0f) + qMat.transpose() * (p2 - p1) - parallelogramOrigin;

	Vec2 opt = parallelogramContainsPt(parallelogramOrigin, parallelogramSideA, parallelogramSideB, Vec2(0.0f, 0.0f));
	if (opt_out)
	{
		*opt_out = opt;
	}

	// sqDist
	return opt.dot(opt) + (p2 - p1).dot(p2 - p1) - (p2 - p1).dot(qMat * qMat.transpose() * (p2 - p1));
}

int manifoldQR(Contact* contacts, const Capsule& bodyA, const Capsule& bodyB)
{
	Vec2 p1 = bodyA.translation() - 0.5f * bodyA.height * bodyA.rotation().col(1); // seg0p0
	Vec2 u1 = bodyA.translation() + 0.5f * bodyA.height * bodyA.rotation().col(1); // seg0p1

	Vec2 p2 = bodyB.translation() - 0.5f * bodyB.height * bodyB.rotation().col(1); // seg1p0
	Vec2 u2 = bodyB.translation() + 0.5f * bodyB.height * bodyB.rotation().col(1); // seg1p1

	Vec2 s1 = u1 - p1;
	Vec2 s2 = u2 - p2;

	float r1 = bodyA.rad;
	float r2 = bodyB.rad;

	Vec2 p1u1[2] = { p1, u1 };
	Vec2 p2u2[2] = { p2, u2 };

	Mat22 qMat, rMat;
	Vec2 opt;

	float sqDist = distanceQR(p1u1, p2u2, &qMat, &rMat, &opt);

	if ((bodyA.rad + bodyB.rad) * (bodyA.rad + bodyB.rad) < sqDist)
	{
		return 0;
	}

	Vec2 rx = opt - qMat.transpose() * (p2 - p1);

	float parallelTest = std::abs(std::abs(s1.dot(s2)) - s1.norm() * s2.norm());

	// line vs line nearest in point - not parallel
	if (std::abs(rMat.col(0).x()) > 1e-4 && std::abs(rMat.col(1).y()) > 1e-4 && parallelTest > 1e-4)
	{
		float x1 = rx.y() / rMat.col(1).y();
		float x2 = (rMat.col(1).y() * rx.x() - rMat.col(1).x() * rx.y()) / (rMat.col(0).x() * rMat.col(1).y());

		if (x2 < 0.0f)
			x2 = 0.0f;

		FeaturePair fp;
		if (x1 <= 0.0f + 1e-4)
		{
			fp.e.inFeature = NEG_Y;
		}
		else if (x1 >= 1.0f - 1e-4)
		{
			fp.e.inFeature = POS_Y;
		}
		else
		{
			fp.e.inFeature = BODY;
		}
		if (x2 <= 0.0f + 1e-4)
		{
			fp.e.outFeature = NEG_Y;
		}
		else if (x2 >= 1.0f - 1e-4)
		{
			fp.e.outFeature = POS_Y;
		}
		else
		{
			fp.e.outFeature = BODY;
		}
		Vec2 nContact = ((p2 + x2 * s2) - (p1 + x1 * s1)).normalized();
		contacts[0].separation = sqrtf(sqDist) - r1 - r2;
		contacts[0].normal = nContact;
		contacts[0].position = 0.5f * (((p1 + x1 * s1) + r1 * nContact) + ((p2 + x2 * s2) - r2 * nContact));
		contacts[0].feature = fp;
		return 1;
	}

	float s1LenSqd = s1.dot(s1);
	float s2LenSqd = s2.dot(s2);

	// point vs point
	if (s1LenSqd < 1e-4 && s2LenSqd < 1e-4)
	{
		FeaturePair fp;
		fp.e.inFeature = BODY;
		fp.e.outFeature = BODY;
		Vec2 nContact = (p2 - p1);
		float nContactNorm = nContact.norm();
		if (nContactNorm > 1e-4)
			nContact /= nContactNorm;
		else
			nContact = Vec2(0.0f, 0.0f);
		contacts[0].separation = sqrtf(sqDist) - r1 - r2;
		contacts[0].normal = nContact;
		contacts[0].position = 0.5f * ((p1 + r1 * nContact) + (p2 - r2 * nContact));
		contacts[0].feature = fp;
		return 1;
	}

	// line vs line - parallel lines
	if (s1LenSqd * s2LenSqd > 1e-4)
	{
		float x1[2];
		float x2[2];
		FeaturePair fp;
		Vec2 nContact;

		if (s1LenSqd > s2LenSqd)
		{
			x1[0] = std::min(std::max(((1.0f / sqrtf(s1LenSqd)) * s1).dot(p2 - p1) / sqrtf(s1LenSqd), 0.0f), 1.0f);
			x1[1] = std::min(std::max(((1.0f / sqrtf(s1LenSqd)) * s1).dot(u2 - p1) / sqrtf(s1LenSqd), 0.0f), 1.0f);

			x2[0] = std::min(std::max(((1.0f / sqrtf(s2LenSqd)) * s2).dot(p1 - p2 + x1[0] * s1) / sqrtf(s2LenSqd), 0.0f), 1.0f);
			x2[1] = std::min(std::max(((1.0f / sqrtf(s2LenSqd)) * s2).dot(p1 - p2 + x1[1] * s1) / sqrtf(s2LenSqd), 0.0f), 1.0f);

		}
		else
		{
			x2[0] = std::min(std::max(((1.0f / sqrtf(s2LenSqd)) * s2).dot(p1 - p2) / sqrtf(s2LenSqd), 0.0f), 1.0f);
			x2[1] = std::min(std::max(((1.0f / sqrtf(s2LenSqd)) * s2).dot(u1 - p2) / sqrtf(s2LenSqd), 0.0f), 1.0f);

			x1[0] = std::min(std::max(((1.0f / sqrtf(s1LenSqd)) * s1).dot(p2 - p1 + x2[0] * s2) / sqrtf(s1LenSqd), 0.0f), 1.0f);
			x1[1] = std::min(std::max(((1.0f / sqrtf(s1LenSqd)) * s1).dot(p2 - p1 + x2[1] * s2) / sqrtf(s1LenSqd), 0.0f), 1.0f);
		}

		nContact = (p2 + 0.5f * (x2[0] + x2[1]) * s2) - (p1 + 0.5f * (x1[0] + x1[1]) * s1);
		float nContactNorm = nContact.norm();
		if (nContactNorm > 1e-4)
			nContact /= nContactNorm;
		else
			nContact = Vec2(0.0f, 0.0f);
		for (int i = 0; i < 2; ++i)
		{
			contacts[i].separation = sqrtf(sqDist) - r1 - r2;
			contacts[i].normal = nContact;
			contacts[i].position = (p1 + x1[i] * s1) + 0.5f * (r1 - r2 + sqrtf(sqDist)) * nContact;
			fp.e.inFeature = BODY;
			fp.e.outFeature = BODY;
			contacts[i].feature = fp;
		}
		if (((x1[0] <= 0.0f + 1e-4 && x1[1] <= 0.0f + 1e-4) || (x1[0] >= 1.0f - 1e-4 && x1[1] >= 1.0f - 1e-4)) && ((x2[0] <= 0.0f + 1e-4 && x2[1] <= 0.0f + 1e-4) || (x2[0] >= 1.0f - 1e-4 && x2[1] >= 1.0f - 1e-4)))
		{
			fp.e.inFeature = BODY;
			fp.e.outFeature = BODY;
			contacts[0].feature = fp;
			return 1;
		}
		else
			return 2;
	}

	// line vs point
	{
		float x1[2];
		float x2[2];
		FeaturePair fp;
		Vec2 nContact;

		if (s1LenSqd > s2LenSqd)
		{
			x1[1] = std::min(std::max(((1.0f / sqrtf(s1LenSqd)) * s1).dot(u2 - p1) / sqrtf(s1LenSqd), 0.0f), 1.0f);
			x1[0] = std::min(std::max(((1.0f / sqrtf(s1LenSqd)) * s1).dot(p2 - p1) / sqrtf(s1LenSqd), 0.0f), 1.0f);

			x2[0] = 0.0f;
			x2[1] = 1.0f;

			if (x1[0] == 0.0f)
			{
				fp.e.inFeature = NEG_Y;
			}
			else if (x1[0] == 1.0f)
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
			x1[0] = 0.0f;
			x1[1] = 1.0f;

			x2[0] = std::min(std::max(((1.0f / sqrtf(s2LenSqd)) * s2).dot(p1 - p2) / sqrtf(s2LenSqd), 0.0f), 1.0f);
			x2[1] = std::min(std::max(((1.0f / sqrtf(s2LenSqd)) * s2).dot(u1 - p2) / sqrtf(s2LenSqd), 0.0f), 1.0f);

			if (x2[0] == 0.0f)
			{
				fp.e.outFeature = NEG_Y;
			}
			else if (x2[0] == 1.0f)
			{
				fp.e.outFeature = POS_Y;
			}
			else
			{
				fp.e.outFeature = BODY;
			}
			fp.e.inFeature = BODY;
		}

		nContact = ((p2 + Vec2(x2[0] * s2.x(), x2[1] * s2.y())) - (p1 + Vec2(x1[0] * s1.x(), x1[1] * s1.y())));
		float nContactNorm = nContact.norm();
		if (nContactNorm > 1e-4)
			nContact /= nContactNorm;
		else
			nContact = Vec2(0.0f, 0.0f);
		contacts[0].separation = sqrtf(sqDist) - r1 - r2;
		contacts[0].normal = nContact;
		contacts[0].position = (p1 + Vec2(x1[0] * s1.x(), x1[1] * s1.y())) + 0.5f * (r1 - r2 + sqrtf(sqDist)) * nContact;
		contacts[0].feature = fp;
		return 1;
	}
}

// 2D rotation matrix around the origin
Mat22 rotation2D(float angle)
{
	Mat22 rotMat = Mat22();
	rotMat(0, 0) = std::cosf(angle);
	rotMat(0, 1) = -std::sinf(angle);
	rotMat(1, 0) = std::sinf(angle);
	rotMat(1, 1) = std::cosf(angle);

	return rotMat;
}

// 2D homogeneous matrix
Mat33 homogeneous2D(Mat22 linear, Vec2 displacement)
{
	Mat33 hMat = Mat33::Identity();
	hMat.block<2, 2>(0, 0) = linear;
	hMat.block<2, 1>(0, 2) = displacement;

	return hMat;
}
