#include <iostream>
#include "include/capsule_contact_manifold_qr.h"

int main()
{
	{
		// line vs line
		Vec2 lineSegA[2] = { Vec2(0.0f, 0.0f), Vec2(-1.0f, 1.0f) };
		Vec2 lineSegB[2] = { Vec2(2.0f, 1.0f), Vec2(-1.0f, 2.0f) };

		float sqDist = distanceQR(lineSegA, lineSegB);
	}

	{
		// line vs line (parallel)
		Vec2 lineSegA[2] = { Vec2(-0.5f, 0.0f), Vec2(0.5f, 0.0f) };
		Vec2 lineSegB[2] = { Vec2(0.25f, 1.0f), Vec2(1.0f, 1.0f) };
		
		float sqDist = distanceQR(lineSegA, lineSegB);
	}

	{
		// line vs point
		Vec2 lineSegA[2] = { Vec2(0.0f, 0.0f), Vec2(2.0f, 2.0f) };
		Vec2 lineSegB[2] = { Vec2(3.0f, 0.0f), Vec2(3.0f, 0.0f) };
		
		float sqDist = distanceQR(lineSegA, lineSegB);
	}

	{
		// point vs point
		Vec2 lineSegA[2] = { Vec2(0.0f, 0.0f), Vec2(0.0f, 0.0f) };
		Vec2 lineSegB[2] = { Vec2(1.0f, 1.0f), Vec2(1.0f, 1.0f) };
		
		float sqDist = distanceQR(lineSegA, lineSegB);
	}
	
	Contact contacts[2];
	int numContacts = 0;

	{
		Mat33 poseA = homogeneous2D(rotation2D(-static_cast<float>(M_PI) / 4.0f), Vec2(1.0f,1.0f));
		float heightA = std::sqrtf(8.0f);
		float radA = 0.75f;
		
		Mat33 poseB = homogeneous2D(Mat22::Identity(), Vec2(3.0f, 0.0f));
		float heightB = 0.0f;
		float radB = 1.5f;

		Capsule bodyA(poseA, heightA, radA);
		Capsule bodyB(poseB, heightB, radB);
		numContacts = manifoldQR(contacts, bodyA, bodyB);
	}

	{
		Mat33 poseA = homogeneous2D(rotation2D(-static_cast<float>(M_PI) / 2.0f), Vec2(1.5f, 0.0f));
		float heightA = 3.0f;
		float radA = 0.75f;
		
		Mat33 poseB = homogeneous2D(rotation2D(-static_cast<float>(M_PI) / 2.0f), Vec2(3.5f, 2.0f));
		float heightB = 5.0f;
		float radB = 1.5f;

		Capsule bodyA(poseA, heightA, radA);
		Capsule bodyB(poseB, heightB, radB);
		numContacts = manifoldQR(contacts, bodyA, bodyB);
	}

	{
		Mat33 poseA = homogeneous2D(rotation2D(-static_cast<float>(M_PI) / 4.0f), Vec2(0.25f, 0.25f));
		float heightA = 1.0f / std::sqrtf(2.0f);
		float radA = 0.5f;

		Mat33 poseB = homogeneous2D(rotation2D(-static_cast<float>(M_PI) / 4.0f), Vec2(1.75f, 1.0f));
		float heightB = std::sqrtf(2.0f);
		float radB = 0.5f;

		Capsule bodyA(poseA, heightA, radA);
		Capsule bodyB(poseB, heightB, radB);
		numContacts = manifoldQR(contacts, bodyA, bodyB);
	}

	std::cout << std::endl;
	return 0;
}
