#include <iostream>

#include <Mathematics/Vector.h>
#include <Mathematics/Segment.h>
#include <Mathematics/Capsule.h>

#include <Mathematics/ArbitraryPrecision.h>
#include <Mathematics/DistSegmentSegment.h>
#include <Mathematics/DistSegmentSegmentExact.h>

#include "include/capsule_contact_manifold_qr.h"

typedef Eigen::Matrix<float, 2, 2> Mat22;
typedef Eigen::Matrix<float, 3, 3> Mat33;
typedef Eigen::Matrix<float, 4, 4> Mat44;
typedef Eigen::Matrix<float, 3, 2> Mat32;
typedef Eigen::Matrix<float, 2, 1> Vec2;
typedef Eigen::Matrix<float, 3, 1> Vec3;

int main()
{
	gte::Vector<3, float> p1 = { 0.0f, 0.0f, 0.0f };
	gte::Vector<3, float> u1 = { -1.0f, 1.0f, 0.0f };
	
	gte::Vector<3, float> p2 = { 2.0f, 1.0f, 0.0f };
	gte::Vector<3, float> u2 = { -1.0f, 2.0f, 0.0f };

	gte::Segment(p1, u1);
	gte::Segment(p2, u2);

	gte::DCPQuery<float, gte::Segment<3, float>, gte::Segment<3, float>> query1;
	auto result1 = query1(p1, u1, p2, u2);
	float distCG = result1.distance;

	Vec3 lineSegA_3D[2] = { Vec3(0.0f, 0.0f, 0.0f), Vec3(-1.0f, 1.0f, 0.0f) };
	Vec3 lineSegB_3D[2] = { Vec3(2.0f, 1.0f, 0.0f), Vec3(-1.0f, 2.0f, 0.0f) };

	float distQR = std::sqrtf(distanceQR(lineSegA_3D, lineSegB_3D));

	return 0;
}
