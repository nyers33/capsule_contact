#include <iostream>
#include "include/capsule_contact_manifold_qr.h"

int main()
{
	{
		// line vs line
		Vec2 lineSegA_2D[2] = { Vec2(0.0f, 0.0f), Vec2(-1.0f, 1.0f) };
		Vec2 lineSegB_2D[2] = { Vec2(2.0f, 1.0f), Vec2(-1.0f, 2.0f) };

		Vec3 lineSegA_3D[2] = { Vec3(0.0f, 0.0f, 0.0f), Vec3(-1.0f, 1.0f, 0.0f) };
		Vec3 lineSegB_3D[2] = { Vec3(2.0f, 1.0f, 0.0f), Vec3(-1.0f, 2.0f, 0.0f) };

		float sqDist_2D = distanceQR<Vec2, Mat22>(lineSegA_2D, lineSegB_2D);
		float sqDist_3D = distanceQR<Vec3, Mat32>(lineSegA_3D, lineSegB_3D);

		std::cout << sqDist_2D << std::endl;
		std::cout << sqDist_3D << std::endl;
	}

	{
		// line vs line (parallel)
		Vec2 lineSegA_2D[2] = { Vec2(-0.5f, 0.0f), Vec2(0.5f, 0.0f) };
		Vec2 lineSegB_2D[2] = { Vec2(0.25f, 1.0f), Vec2(1.0f, 1.0f) };

		Vec3 lineSegA_3D[2] = { Vec3(-0.5f, 0.0f, 0.0f), Vec3(0.5f, 0.0f, 0.0f) };
		Vec3 lineSegB_3D[2] = { Vec3(0.25f, 1.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f) };

		float sqDist_2D = distanceQR<Vec2, Mat22>(lineSegA_2D, lineSegB_2D);
		float sqDist_3D = distanceQR<Vec3, Mat32>(lineSegA_3D, lineSegB_3D);

		std::cout << sqDist_2D << std::endl;
		std::cout << sqDist_3D << std::endl;
	}

	{
		// line vs point
		Vec2 lineSegA_2D[2] = { Vec2(0.0f, 0.0f), Vec2(2.0f, 2.0f) };
		Vec2 lineSegB_2D[2] = { Vec2(3.0f, 0.0f), Vec2(3.0f, 0.0f) };

		Vec3 lineSegA_3D[2] = { Vec3(0.0f, 0.0f, 0.0f), Vec3(2.0f, 2.0f, 0.0f) };
		Vec3 lineSegB_3D[2] = { Vec3(3.0f, 0.0f, 0.0f), Vec3(3.0f, 0.0f, 0.0f) };

		float sqDist_2D = distanceQR<Vec2, Mat22>(lineSegA_2D, lineSegB_2D);
		float sqDist_3D = distanceQR<Vec3, Mat32>(lineSegA_3D, lineSegB_3D);

		std::cout << sqDist_2D << std::endl;
		std::cout << sqDist_3D << std::endl;
	}

	{
		// point vs point
		Vec2 lineSegA_2D[2] = { Vec2(0.0f, 0.0f), Vec2(0.0f, 0.0f) };
		Vec2 lineSegB_2D[2] = { Vec2(1.0f, 1.0f), Vec2(1.0f, 1.0f) };

		Vec3 lineSegA_3D[2] = { Vec3(0.0f, 0.0f, 0.0f), Vec3(0.0f, 0.0f, 0.0f) };
		Vec3 lineSegB_3D[2] = { Vec3(1.0f, 1.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f) };

		float sqDist_2D = distanceQR<Vec2, Mat22>(lineSegA_2D, lineSegB_2D);
		float sqDist_3D = distanceQR<Vec3, Mat32>(lineSegA_3D, lineSegB_3D);

		std::cout << sqDist_2D << std::endl;
		std::cout << sqDist_3D << std::endl;
	}
	
	{
		Contact<Mat44> contacts[2];
		int numContacts = 0;

		auto vectest = rotationAlign3D(Vec3(0.0f, 1.0f, 0.0f), Vec3(-1.0f, 1.0f, 1.0f)) * Vec3(0.0f, 1.0f, 0.0f);
		std::cout << vectest.normalized() << std::endl;
		std::cout << Vec3(-1.0f, 1.0f, 1.0f).normalized() << std::endl;

		auto poseA = homogeneous3D( rotationAlign3D(Vec3(0.0f, 1.0f, 0.0f), Vec3(-1.0f, 1.0f, 1.0f)), Vec3(-0.5f, 0.5f, 0.5f) );
		float heightA = std::sqrtf(3.0f);
		float radA = 1.0f;
		
		auto poseB = homogeneous3D(rotationAlign3D(Vec3(0.0f, 1.0f, 0.0f), Vec3(-4.0f, 1.0f, -1.5f)), Vec3(0.0f, 1.5f, 2.25f));
		float heightB = std::sqrtf(77.0f) / 2.0f;
		float radB = 1.0f;

		Capsule bodyA(poseA, heightA, radA);
		Capsule bodyB(poseB, heightB, radB);

		numContacts = manifoldQR(contacts, bodyA, bodyB);
	}

	{
		Contact<Mat44> contacts[2];
		int numContacts = 0;

		auto vectest = rotationAlign3D(Vec3(0.0f, 1.0f, 0.0f), Vec3(-1.5f, -1.5f, 0.0f)) * Vec3(0.0f, 1.0, 0.0f);
		std::cout << vectest.normalized() << std::endl;
		std::cout << Vec3(-1.5f, -1.5f, 0.0f).normalized() << std::endl;

		auto poseA = homogeneous3D(rotationAlign3D(Vec3(0.0f, 1.0f, 0.0f), Vec3(-1.5f, -1.5f, 0.0f)), Vec3(0.25f, -0.25f, 0.0f));
		float heightA = 3.0f / std::sqrtf(2.0f);
		float radA = 1.0f;

		auto poseB = homogeneous3D(rotationAlign3D(Vec3(0.0f, 1.0f, 0.0f), Vec3(-1.0f, 2.0f, 0.0f)), Vec3(0.0f, 0.0f, 0.0f));
		float heightB = std::sqrtf(5.0f);
		float radB = 1.0f;

		Capsule bodyA(poseA, heightA, radA);
		Capsule bodyB(poseB, heightB, radB);

		numContacts = manifoldQR(contacts, bodyA, bodyB);
	}

	{
		Contact<Mat44> contacts[2];
		int numContacts = 0;

		auto vectest = rotationAlign3D(Vec3(0.0f, 1.0f, 0.0f), Vec3(-2.0f, -2.0f, 0.0f)) * Vec3(0.0f, 1.0, 0.0f);
		std::cout << vectest.normalized() << std::endl;
		std::cout << Vec3(-2.0f, -2.0f, 0.0f).normalized() << std::endl;

		auto poseA = homogeneous3D(rotationAlign3D(Vec3(0.0f, 1.0f, 0.0f), Vec3(-2.0f, -2.0f, 0.0f)), Vec3(-1.0f, 0.0f, 0.0f));
		float heightA = 2.0f * std::sqrtf(2.0f);
		float radA = 1.0f;

		auto poseB = homogeneous3D(rotationAlign3D(Vec3(0.0f, 1.0f, 0.0f), Vec3(0.0f, 0.0f, -1.0f)), Vec3(0.0f, 0.0f, 0.0f));
		float heightB = 1.0f;
		float radB = 1.0f;

		Capsule bodyA(poseA, heightA, radA);
		Capsule bodyB(poseB, heightB, radB);

		numContacts = manifoldQR(contacts, bodyA, bodyB);
	}

	std::cout << std::endl;

	return 0;
}
