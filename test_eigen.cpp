#include <iostream>
#include "include/capsule_contact_manifold_qr.h"

typedef Eigen::Matrix<float, 2, 2> Mat22;
typedef Eigen::Matrix<float, 3, 3> Mat33;
typedef Eigen::Matrix<float, 4, 4> Mat44;
typedef Eigen::Matrix<float, 3, 2> Mat32;
typedef Eigen::Matrix<float, 2, 1> Vec2;
typedef Eigen::Matrix<float, 3, 1> Vec3;

int main()
{
	{
		// line vs line
		Vec2 lineSegA_2D[2] = { Vec2(0.0f, 0.0f), Vec2(-1.0f, 1.0f) };
		Vec2 lineSegB_2D[2] = { Vec2(2.0f, 1.0f), Vec2(-1.0f, 2.0f) };

		Vec3 lineSegA_3D[2] = { Vec3(0.0f, 0.0f, 0.0f), Vec3(-1.0f, 1.0f, 0.0f) };
		Vec3 lineSegB_3D[2] = { Vec3(2.0f, 1.0f, 0.0f), Vec3(-1.0f, 2.0f, 0.0f) };

		float sqDist_2D = distanceQR(lineSegA_2D, lineSegB_2D);
		float sqDist_3D = distanceQR(lineSegA_3D, lineSegB_3D);

		std::cout << sqDist_2D << std::endl;
		std::cout << sqDist_3D << std::endl;
	}

	{
		// line vs line (parallel)
		Vec2 lineSegA_2D[2] = { Vec2(-0.5f, 0.0f), Vec2(0.5f, 0.0f) };
		Vec2 lineSegB_2D[2] = { Vec2(0.25f, 1.0f), Vec2(1.0f, 1.0f) };

		Vec3 lineSegA_3D[2] = { Vec3(-0.5f, 0.0f, 0.0f), Vec3(0.5f, 0.0f, 0.0f) };
		Vec3 lineSegB_3D[2] = { Vec3(0.25f, 1.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f) };

		float sqDist_2D = distanceQR(lineSegA_2D, lineSegB_2D);
		float sqDist_3D = distanceQR(lineSegA_3D, lineSegB_3D);

		std::cout << sqDist_2D << std::endl;
		std::cout << sqDist_3D << std::endl;
	}

	{
		// line vs point
		Vec2 lineSegA_2D[2] = { Vec2(0.0f, 0.0f), Vec2(2.0f, 2.0f) };
		Vec2 lineSegB_2D[2] = { Vec2(3.0f, 0.0f), Vec2(3.0f, 0.0f) };

		Vec3 lineSegA_3D[2] = { Vec3(0.0f, 0.0f, 0.0f), Vec3(2.0f, 2.0f, 0.0f) };
		Vec3 lineSegB_3D[2] = { Vec3(3.0f, 0.0f, 0.0f), Vec3(3.0f, 0.0f, 0.0f) };

		float sqDist_2D = distanceQR(lineSegA_2D, lineSegB_2D);
		float sqDist_3D = distanceQR(lineSegA_3D, lineSegB_3D);

		std::cout << sqDist_2D << std::endl;
		std::cout << sqDist_3D << std::endl;
	}

	{
		// point vs point
		Vec2 lineSegA_2D[2] = { Vec2(0.0f, 0.0f), Vec2(0.0f, 0.0f) };
		Vec2 lineSegB_2D[2] = { Vec2(1.0f, 1.0f), Vec2(1.0f, 1.0f) };

		Vec3 lineSegA_3D[2] = { Vec3(0.0f, 0.0f, 0.0f), Vec3(0.0f, 0.0f, 0.0f) };
		Vec3 lineSegB_3D[2] = { Vec3(1.0f, 1.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f) };

		float sqDist_2D = distanceQR(lineSegA_2D, lineSegB_2D);
		float sqDist_3D = distanceQR(lineSegA_3D, lineSegB_3D);

		std::cout << sqDist_2D << std::endl;
		std::cout << sqDist_3D << std::endl;
	}
	
	{
		Contact<float, 3> contacts[2];
		int numContacts = 0;

		Mat44 poseA = homogeneous3D( rotationAlign3D(Vec3(0.0f, 1.0f, 0.0f), Vec3(-1.0f, 1.0f, -1.0f)), Vec3(-0.5f, 0.5f, -0.5f) );
		float heightA = std::sqrtf(3.0f);
		float radA = 2.0f;
		
		Mat44 poseB = homogeneous3D(rotationAlign3D(Vec3(0.0f, 1.0f, 0.0f), Vec3(-3.0f, 2.0f, -2.0f)), Vec3(0.5f, 2.0f, 2.0f));
		float heightB = std::sqrtf(17.0f);
		float radB = 2.0f;

		Capsule<float, 3> bodyA(poseA, heightA, radA);
		Capsule<float, 3> bodyB(poseB, heightB, radB);

		numContacts = manifoldQR(contacts, bodyA, bodyB);
	}

	{
		Contact<float, 3> contacts[2];
		int numContacts = 0;

		Mat44 poseA = homogeneous3D(rotationAlign3D(Vec3(0.0f, 1.0f, 0.0f), Vec3(-1.5f, -1.5f, 0.0f)), Vec3(0.25f, -0.25f, 0.0f));
		float heightA = 3.0f / std::sqrtf(2.0f);
		float radA = 1.0f;

		Mat44 poseB = homogeneous3D(rotationAlign3D(Vec3(0.0f, 1.0f, 0.0f), Vec3(-1.0f, 2.0f, 0.0f)), Vec3(0.0f, 0.0f, 0.0f));
		float heightB = std::sqrtf(5.0f);
		float radB = 1.0f;

		Capsule<float, 3> bodyA(poseA, heightA, radA);
		Capsule<float, 3> bodyB(poseB, heightB, radB);

		numContacts = manifoldQR(contacts, bodyA, bodyB);
	}

	{
		Contact<float, 3> contacts[2];
		int numContacts = 0;

		Mat44 poseA = homogeneous3D(rotationAlign3D(Vec3(0.0f, 1.0f, 0.0f), Vec3(-2.0f, -2.0f, 0.0f)), Vec3(-1.0f, 0.0f, 0.0f));
		float heightA = 2.0f * std::sqrtf(2.0f);
		float radA = 1.0f;

		Mat44 poseB = homogeneous3D(rotationAlign3D(Vec3(0.0f, 1.0f, 0.0f), Vec3(0.0f, 0.0f, -1.0f)), Vec3(0.0f, 0.0f, 0.0f));
		float heightB = 1.0f;
		float radB = 1.0f;

		Capsule<float, 3> bodyA(poseA, heightA, radA);
		Capsule<float, 3> bodyB(poseB, heightB, radB);

		numContacts = manifoldQR(contacts, bodyA, bodyB);
	}
	
	std::cout << std::endl;

	return 0;
}
