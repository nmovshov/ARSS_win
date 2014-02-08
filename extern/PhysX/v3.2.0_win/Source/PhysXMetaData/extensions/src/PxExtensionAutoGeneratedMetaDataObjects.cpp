// This code contains NVIDIA Confidential Information and is disclosed to you 
// under a form of NVIDIA software license agreement provided separately to you.
//
// Notice
// NVIDIA Corporation and its licensors retain all intellectual property and
// proprietary rights in and to this software and related documentation and 
// any modifications thereto. Any use, reproduction, disclosure, or 
// distribution of this software and related documentation without an express 
// license agreement from NVIDIA Corporation is strictly prohibited.
// 
// ALL NVIDIA DESIGN SPECIFICATIONS, CODE ARE PROVIDED "AS IS.". NVIDIA MAKES
// NO WARRANTIES, EXPRESSED, IMPLIED, STATUTORY, OR OTHERWISE WITH RESPECT TO
// THE MATERIALS, AND EXPRESSLY DISCLAIMS ALL IMPLIED WARRANTIES OF NONINFRINGEMENT,
// MERCHANTABILITY, AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Information and code furnished is believed to be accurate and reliable.
// However, NVIDIA Corporation assumes no responsibility for the consequences of use of such
// information or for any infringement of patents or other rights of third parties that may
// result from its use. No license is granted by implication or otherwise under any patent
// or patent rights of NVIDIA Corporation. Details are subject to change without notice.
// This code supersedes and replaces all information previously supplied.
// NVIDIA Corporation products are not authorized for use as critical
// components in life support devices or systems without express written approval of
// NVIDIA Corporation.
//
// Copyright (c) 2008-2012 NVIDIA Corporation. All rights reserved.
// Copyright (c) 2004-2008 AGEIA Technologies, Inc. All rights reserved.
// Copyright (c) 2001-2004 NovodeX AG. All rights reserved.

// This code is auto-generated by the PhysX Clang metadata generator.  Do not edit or be
// prepared for your edits to be quietly ignored next time the clang metadata generator is
// run.  You can find the most recent version of clang metadata generator by contacting
// Chris Nuernberger <chrisn@nvidia.com> or Dilip or Adam.
// The source code for the generate was at one time checked into:
// physx/PhysXMetaDataGenerator/llvm/tools/clang/lib/Frontend/PhysXMetaDataAction.cpp
#include "PxExtensionMetaDataObjects.h"
#include "PxPhysicsAPI.h"
#include "PsIntrinsics.h"
#include "PxMetaDataCppPrefix.h"
using namespace physx;
#include "PxExtensionsAPI.h"
void setPxJoint_Actors( PxJoint* inObj, PxRigidActor * inArg0, PxRigidActor * inArg1 ) { inObj->setActors( inArg0, inArg1 ); }
void getPxJoint_Actors( const PxJoint* inObj, PxRigidActor *& inArg0, PxRigidActor *& inArg1 ) { inObj->getActors( inArg0, inArg1 ); }
void setPxJoint_LocalPose( PxJoint* inObj, PxJointActorIndex::Enum inIndex, PxTransform inArg ){ inObj->setLocalPose( inIndex, inArg ); }
PxTransform getPxJoint_LocalPose( const PxJoint* inObj, PxJointActorIndex::Enum inIndex ) { return inObj->getLocalPose( inIndex ); }
void setPxJoint_BreakForce( PxJoint* inObj, PxReal inArg0, PxReal inArg1 ) { inObj->setBreakForce( inArg0, inArg1 ); }
void getPxJoint_BreakForce( const PxJoint* inObj, PxReal& inArg0, PxReal& inArg1 ) { inObj->getBreakForce( inArg0, inArg1 ); }
void setPxJoint_ConstraintFlags( PxJoint* inObj, PxConstraintFlags inArg){ inObj->setConstraintFlags( inArg ); }
PxConstraintFlags getPxJoint_ConstraintFlags( const PxJoint* inObj ) { return inObj->getConstraintFlags(); }
PxConstraint * getPxJoint_Constraint( const PxJoint* inObj ) { return inObj->getConstraint(); }
void setPxJoint_Name( PxJoint* inObj, const char * inArg){ inObj->setName( inArg ); }
const char * getPxJoint_Name( const PxJoint* inObj ) { return inObj->getName(); }
PxScene * getPxJoint_Scene( const PxJoint* inObj ) { return inObj->getScene(); }
PxJointType::Enum getPxJoint_Type( const PxJoint* inObj ) { return inObj->getType(); }
inline void * getPxJointUserData( const PxJoint* inOwner ) { return inOwner->userData; }
inline void setPxJointUserData( PxJoint* inOwner, void * inData) { inOwner->userData = inData; }
 PxJointGeneratedInfo::PxJointGeneratedInfo()
	: Actors( "Actors", "actor0", "actor1", setPxJoint_Actors, getPxJoint_Actors)
	, LocalPose( "LocalPose", setPxJoint_LocalPose, getPxJoint_LocalPose)
	, BreakForce( "BreakForce", "force", "torque", setPxJoint_BreakForce, getPxJoint_BreakForce)
	, ConstraintFlags( "ConstraintFlags", setPxJoint_ConstraintFlags, getPxJoint_ConstraintFlags)
	, Constraint( "Constraint", getPxJoint_Constraint)
	, Name( "Name", setPxJoint_Name, getPxJoint_Name)
	, Scene( "Scene", getPxJoint_Scene)
	, Type( "Type", getPxJoint_Type)
	, UserData( "UserData", setPxJointUserData, getPxJointUserData )
{}
 PxJointGeneratedValues::PxJointGeneratedValues( const PxJoint* inSource )
		:ConstraintFlags( getPxJoint_ConstraintFlags( inSource ) )
		,Constraint( getPxJoint_Constraint( inSource ) )
		,Name( getPxJoint_Name( inSource ) )
		,Scene( getPxJoint_Scene( inSource ) )
		,Type( getPxJoint_Type( inSource ) )
		,UserData( inSource->userData )
{
	getPxJoint_Actors( inSource, Actors[0], Actors[1] );
		for ( PxU32 idx = 0; idx < static_cast<PxU32>( physx::PxJointActorIndex::COUNT ); ++idx )
		LocalPose[idx] = getPxJoint_LocalPose( inSource, static_cast< PxJointActorIndex::Enum >( idx ) );
	getPxJoint_BreakForce( inSource, BreakForce[0], BreakForce[1] );
}
void setPxFixedJoint_ProjectionLinearTolerance( PxFixedJoint* inObj, PxReal inArg){ inObj->setProjectionLinearTolerance( inArg ); }
PxReal getPxFixedJoint_ProjectionLinearTolerance( const PxFixedJoint* inObj ) { return inObj->getProjectionLinearTolerance(); }
void setPxFixedJoint_ProjectionAngularTolerance( PxFixedJoint* inObj, PxReal inArg){ inObj->setProjectionAngularTolerance( inArg ); }
PxReal getPxFixedJoint_ProjectionAngularTolerance( const PxFixedJoint* inObj ) { return inObj->getProjectionAngularTolerance(); }
const char * getPxFixedJoint_ConcreteTypeName( const PxFixedJoint* inObj ) { return inObj->getConcreteTypeName(); }
 PxFixedJointGeneratedInfo::PxFixedJointGeneratedInfo()
	: ProjectionLinearTolerance( "ProjectionLinearTolerance", setPxFixedJoint_ProjectionLinearTolerance, getPxFixedJoint_ProjectionLinearTolerance)
	, ProjectionAngularTolerance( "ProjectionAngularTolerance", setPxFixedJoint_ProjectionAngularTolerance, getPxFixedJoint_ProjectionAngularTolerance)
	, ConcreteTypeName( "ConcreteTypeName", getPxFixedJoint_ConcreteTypeName)
{}
 PxFixedJointGeneratedValues::PxFixedJointGeneratedValues( const PxFixedJoint* inSource )
		:PxJointGeneratedValues( inSource )
		,ProjectionLinearTolerance( getPxFixedJoint_ProjectionLinearTolerance( inSource ) )
		,ProjectionAngularTolerance( getPxFixedJoint_ProjectionAngularTolerance( inSource ) )
		,ConcreteTypeName( getPxFixedJoint_ConcreteTypeName( inSource ) )
{
}
void setPxRevoluteJoint_Limit( PxRevoluteJoint* inObj, const PxJointLimitPair & inArg){ inObj->setLimit( inArg ); }
PxJointLimitPair getPxRevoluteJoint_Limit( const PxRevoluteJoint* inObj ) { return inObj->getLimit(); }
void setPxRevoluteJoint_DriveVelocity( PxRevoluteJoint* inObj, PxReal inArg){ inObj->setDriveVelocity( inArg ); }
PxReal getPxRevoluteJoint_DriveVelocity( const PxRevoluteJoint* inObj ) { return inObj->getDriveVelocity(); }
void setPxRevoluteJoint_DriveForceLimit( PxRevoluteJoint* inObj, PxReal inArg){ inObj->setDriveForceLimit( inArg ); }
PxReal getPxRevoluteJoint_DriveForceLimit( const PxRevoluteJoint* inObj ) { return inObj->getDriveForceLimit(); }
void setPxRevoluteJoint_DriveGearRatio( PxRevoluteJoint* inObj, PxReal inArg){ inObj->setDriveGearRatio( inArg ); }
PxReal getPxRevoluteJoint_DriveGearRatio( const PxRevoluteJoint* inObj ) { return inObj->getDriveGearRatio(); }
void setPxRevoluteJoint_RevoluteJointFlags( PxRevoluteJoint* inObj, PxRevoluteJointFlags inArg){ inObj->setRevoluteJointFlags( inArg ); }
PxRevoluteJointFlags getPxRevoluteJoint_RevoluteJointFlags( const PxRevoluteJoint* inObj ) { return inObj->getRevoluteJointFlags(); }
void setPxRevoluteJoint_ProjectionLinearTolerance( PxRevoluteJoint* inObj, PxReal inArg){ inObj->setProjectionLinearTolerance( inArg ); }
PxReal getPxRevoluteJoint_ProjectionLinearTolerance( const PxRevoluteJoint* inObj ) { return inObj->getProjectionLinearTolerance(); }
void setPxRevoluteJoint_ProjectionAngularTolerance( PxRevoluteJoint* inObj, PxReal inArg){ inObj->setProjectionAngularTolerance( inArg ); }
PxReal getPxRevoluteJoint_ProjectionAngularTolerance( const PxRevoluteJoint* inObj ) { return inObj->getProjectionAngularTolerance(); }
const char * getPxRevoluteJoint_ConcreteTypeName( const PxRevoluteJoint* inObj ) { return inObj->getConcreteTypeName(); }
 PxRevoluteJointGeneratedInfo::PxRevoluteJointGeneratedInfo()
	: Limit( "Limit", setPxRevoluteJoint_Limit, getPxRevoluteJoint_Limit)
	, DriveVelocity( "DriveVelocity", setPxRevoluteJoint_DriveVelocity, getPxRevoluteJoint_DriveVelocity)
	, DriveForceLimit( "DriveForceLimit", setPxRevoluteJoint_DriveForceLimit, getPxRevoluteJoint_DriveForceLimit)
	, DriveGearRatio( "DriveGearRatio", setPxRevoluteJoint_DriveGearRatio, getPxRevoluteJoint_DriveGearRatio)
	, RevoluteJointFlags( "RevoluteJointFlags", setPxRevoluteJoint_RevoluteJointFlags, getPxRevoluteJoint_RevoluteJointFlags)
	, ProjectionLinearTolerance( "ProjectionLinearTolerance", setPxRevoluteJoint_ProjectionLinearTolerance, getPxRevoluteJoint_ProjectionLinearTolerance)
	, ProjectionAngularTolerance( "ProjectionAngularTolerance", setPxRevoluteJoint_ProjectionAngularTolerance, getPxRevoluteJoint_ProjectionAngularTolerance)
	, ConcreteTypeName( "ConcreteTypeName", getPxRevoluteJoint_ConcreteTypeName)
{}
 PxRevoluteJointGeneratedValues::PxRevoluteJointGeneratedValues( const PxRevoluteJoint* inSource )
		:PxJointGeneratedValues( inSource )
		,Limit( getPxRevoluteJoint_Limit( inSource ) )
		,DriveVelocity( getPxRevoluteJoint_DriveVelocity( inSource ) )
		,DriveForceLimit( getPxRevoluteJoint_DriveForceLimit( inSource ) )
		,DriveGearRatio( getPxRevoluteJoint_DriveGearRatio( inSource ) )
		,RevoluteJointFlags( getPxRevoluteJoint_RevoluteJointFlags( inSource ) )
		,ProjectionLinearTolerance( getPxRevoluteJoint_ProjectionLinearTolerance( inSource ) )
		,ProjectionAngularTolerance( getPxRevoluteJoint_ProjectionAngularTolerance( inSource ) )
		,ConcreteTypeName( getPxRevoluteJoint_ConcreteTypeName( inSource ) )
{
}
void setPxPrismaticJoint_Limit( PxPrismaticJoint* inObj, const PxJointLimitPair & inArg){ inObj->setLimit( inArg ); }
PxJointLimitPair getPxPrismaticJoint_Limit( const PxPrismaticJoint* inObj ) { return inObj->getLimit(); }
void setPxPrismaticJoint_PrismaticJointFlags( PxPrismaticJoint* inObj, PxPrismaticJointFlags inArg){ inObj->setPrismaticJointFlags( inArg ); }
PxPrismaticJointFlags getPxPrismaticJoint_PrismaticJointFlags( const PxPrismaticJoint* inObj ) { return inObj->getPrismaticJointFlags(); }
void setPxPrismaticJoint_ProjectionLinearTolerance( PxPrismaticJoint* inObj, PxReal inArg){ inObj->setProjectionLinearTolerance( inArg ); }
PxReal getPxPrismaticJoint_ProjectionLinearTolerance( const PxPrismaticJoint* inObj ) { return inObj->getProjectionLinearTolerance(); }
void setPxPrismaticJoint_ProjectionAngularTolerance( PxPrismaticJoint* inObj, PxReal inArg){ inObj->setProjectionAngularTolerance( inArg ); }
PxReal getPxPrismaticJoint_ProjectionAngularTolerance( const PxPrismaticJoint* inObj ) { return inObj->getProjectionAngularTolerance(); }
const char * getPxPrismaticJoint_ConcreteTypeName( const PxPrismaticJoint* inObj ) { return inObj->getConcreteTypeName(); }
 PxPrismaticJointGeneratedInfo::PxPrismaticJointGeneratedInfo()
	: Limit( "Limit", setPxPrismaticJoint_Limit, getPxPrismaticJoint_Limit)
	, PrismaticJointFlags( "PrismaticJointFlags", setPxPrismaticJoint_PrismaticJointFlags, getPxPrismaticJoint_PrismaticJointFlags)
	, ProjectionLinearTolerance( "ProjectionLinearTolerance", setPxPrismaticJoint_ProjectionLinearTolerance, getPxPrismaticJoint_ProjectionLinearTolerance)
	, ProjectionAngularTolerance( "ProjectionAngularTolerance", setPxPrismaticJoint_ProjectionAngularTolerance, getPxPrismaticJoint_ProjectionAngularTolerance)
	, ConcreteTypeName( "ConcreteTypeName", getPxPrismaticJoint_ConcreteTypeName)
{}
 PxPrismaticJointGeneratedValues::PxPrismaticJointGeneratedValues( const PxPrismaticJoint* inSource )
		:PxJointGeneratedValues( inSource )
		,Limit( getPxPrismaticJoint_Limit( inSource ) )
		,PrismaticJointFlags( getPxPrismaticJoint_PrismaticJointFlags( inSource ) )
		,ProjectionLinearTolerance( getPxPrismaticJoint_ProjectionLinearTolerance( inSource ) )
		,ProjectionAngularTolerance( getPxPrismaticJoint_ProjectionAngularTolerance( inSource ) )
		,ConcreteTypeName( getPxPrismaticJoint_ConcreteTypeName( inSource ) )
{
}
void setPxSphericalJoint_LimitCone( PxSphericalJoint* inObj, const PxJointLimitCone & inArg){ inObj->setLimitCone( inArg ); }
PxJointLimitCone getPxSphericalJoint_LimitCone( const PxSphericalJoint* inObj ) { return inObj->getLimitCone(); }
void setPxSphericalJoint_SphericalJointFlags( PxSphericalJoint* inObj, PxSphericalJointFlags inArg){ inObj->setSphericalJointFlags( inArg ); }
PxSphericalJointFlags getPxSphericalJoint_SphericalJointFlags( const PxSphericalJoint* inObj ) { return inObj->getSphericalJointFlags(); }
void setPxSphericalJoint_ProjectionLinearTolerance( PxSphericalJoint* inObj, PxReal inArg){ inObj->setProjectionLinearTolerance( inArg ); }
PxReal getPxSphericalJoint_ProjectionLinearTolerance( const PxSphericalJoint* inObj ) { return inObj->getProjectionLinearTolerance(); }
const char * getPxSphericalJoint_ConcreteTypeName( const PxSphericalJoint* inObj ) { return inObj->getConcreteTypeName(); }
 PxSphericalJointGeneratedInfo::PxSphericalJointGeneratedInfo()
	: LimitCone( "LimitCone", setPxSphericalJoint_LimitCone, getPxSphericalJoint_LimitCone)
	, SphericalJointFlags( "SphericalJointFlags", setPxSphericalJoint_SphericalJointFlags, getPxSphericalJoint_SphericalJointFlags)
	, ProjectionLinearTolerance( "ProjectionLinearTolerance", setPxSphericalJoint_ProjectionLinearTolerance, getPxSphericalJoint_ProjectionLinearTolerance)
	, ConcreteTypeName( "ConcreteTypeName", getPxSphericalJoint_ConcreteTypeName)
{}
 PxSphericalJointGeneratedValues::PxSphericalJointGeneratedValues( const PxSphericalJoint* inSource )
		:PxJointGeneratedValues( inSource )
		,LimitCone( getPxSphericalJoint_LimitCone( inSource ) )
		,SphericalJointFlags( getPxSphericalJoint_SphericalJointFlags( inSource ) )
		,ProjectionLinearTolerance( getPxSphericalJoint_ProjectionLinearTolerance( inSource ) )
		,ConcreteTypeName( getPxSphericalJoint_ConcreteTypeName( inSource ) )
{
}
void setPxDistanceJoint_MinDistance( PxDistanceJoint* inObj, PxReal inArg){ inObj->setMinDistance( inArg ); }
PxReal getPxDistanceJoint_MinDistance( const PxDistanceJoint* inObj ) { return inObj->getMinDistance(); }
void setPxDistanceJoint_MaxDistance( PxDistanceJoint* inObj, PxReal inArg){ inObj->setMaxDistance( inArg ); }
PxReal getPxDistanceJoint_MaxDistance( const PxDistanceJoint* inObj ) { return inObj->getMaxDistance(); }
void setPxDistanceJoint_Tolerance( PxDistanceJoint* inObj, PxReal inArg){ inObj->setTolerance( inArg ); }
PxReal getPxDistanceJoint_Tolerance( const PxDistanceJoint* inObj ) { return inObj->getTolerance(); }
void setPxDistanceJoint_Spring( PxDistanceJoint* inObj, PxReal inArg){ inObj->setSpring( inArg ); }
PxReal getPxDistanceJoint_Spring( const PxDistanceJoint* inObj ) { return inObj->getSpring(); }
void setPxDistanceJoint_Damping( PxDistanceJoint* inObj, PxReal inArg){ inObj->setDamping( inArg ); }
PxReal getPxDistanceJoint_Damping( const PxDistanceJoint* inObj ) { return inObj->getDamping(); }
void setPxDistanceJoint_DistanceJointFlags( PxDistanceJoint* inObj, PxDistanceJointFlags inArg){ inObj->setDistanceJointFlags( inArg ); }
PxDistanceJointFlags getPxDistanceJoint_DistanceJointFlags( const PxDistanceJoint* inObj ) { return inObj->getDistanceJointFlags(); }
const char * getPxDistanceJoint_ConcreteTypeName( const PxDistanceJoint* inObj ) { return inObj->getConcreteTypeName(); }
 PxDistanceJointGeneratedInfo::PxDistanceJointGeneratedInfo()
	: MinDistance( "MinDistance", setPxDistanceJoint_MinDistance, getPxDistanceJoint_MinDistance)
	, MaxDistance( "MaxDistance", setPxDistanceJoint_MaxDistance, getPxDistanceJoint_MaxDistance)
	, Tolerance( "Tolerance", setPxDistanceJoint_Tolerance, getPxDistanceJoint_Tolerance)
	, Spring( "Spring", setPxDistanceJoint_Spring, getPxDistanceJoint_Spring)
	, Damping( "Damping", setPxDistanceJoint_Damping, getPxDistanceJoint_Damping)
	, DistanceJointFlags( "DistanceJointFlags", setPxDistanceJoint_DistanceJointFlags, getPxDistanceJoint_DistanceJointFlags)
	, ConcreteTypeName( "ConcreteTypeName", getPxDistanceJoint_ConcreteTypeName)
{}
 PxDistanceJointGeneratedValues::PxDistanceJointGeneratedValues( const PxDistanceJoint* inSource )
		:PxJointGeneratedValues( inSource )
		,MinDistance( getPxDistanceJoint_MinDistance( inSource ) )
		,MaxDistance( getPxDistanceJoint_MaxDistance( inSource ) )
		,Tolerance( getPxDistanceJoint_Tolerance( inSource ) )
		,Spring( getPxDistanceJoint_Spring( inSource ) )
		,Damping( getPxDistanceJoint_Damping( inSource ) )
		,DistanceJointFlags( getPxDistanceJoint_DistanceJointFlags( inSource ) )
		,ConcreteTypeName( getPxDistanceJoint_ConcreteTypeName( inSource ) )
{
}
void setPxD6Joint_Motion( PxD6Joint* inObj, PxD6Axis::Enum inIndex, PxD6Motion::Enum inArg ){ inObj->setMotion( inIndex, inArg ); }
PxD6Motion::Enum getPxD6Joint_Motion( const PxD6Joint* inObj, PxD6Axis::Enum inIndex ) { return inObj->getMotion( inIndex ); }
void setPxD6Joint_LinearLimit( PxD6Joint* inObj, const PxJointLimit & inArg){ inObj->setLinearLimit( inArg ); }
PxJointLimit getPxD6Joint_LinearLimit( const PxD6Joint* inObj ) { return inObj->getLinearLimit(); }
void setPxD6Joint_TwistLimit( PxD6Joint* inObj, const PxJointLimitPair & inArg){ inObj->setTwistLimit( inArg ); }
PxJointLimitPair getPxD6Joint_TwistLimit( const PxD6Joint* inObj ) { return inObj->getTwistLimit(); }
void setPxD6Joint_SwingLimit( PxD6Joint* inObj, const PxJointLimitCone & inArg){ inObj->setSwingLimit( inArg ); }
PxJointLimitCone getPxD6Joint_SwingLimit( const PxD6Joint* inObj ) { return inObj->getSwingLimit(); }
void setPxD6Joint_Drive( PxD6Joint* inObj, PxD6Drive::Enum inIndex, PxD6JointDrive inArg ){ inObj->setDrive( inIndex, inArg ); }
PxD6JointDrive getPxD6Joint_Drive( const PxD6Joint* inObj, PxD6Drive::Enum inIndex ) { return inObj->getDrive( inIndex ); }
void setPxD6Joint_DrivePosition( PxD6Joint* inObj, const PxTransform & inArg){ inObj->setDrivePosition( inArg ); }
PxTransform getPxD6Joint_DrivePosition( const PxD6Joint* inObj ) { return inObj->getDrivePosition(); }
void setPxD6Joint_DriveVelocity( PxD6Joint* inObj, PxVec3 inArg0, PxVec3 inArg1 ) { inObj->setDriveVelocity( inArg0, inArg1 ); }
void getPxD6Joint_DriveVelocity( const PxD6Joint* inObj, PxVec3& inArg0, PxVec3& inArg1 ) { inObj->getDriveVelocity( inArg0, inArg1 ); }
void setPxD6Joint_ProjectionLinearTolerance( PxD6Joint* inObj, PxReal inArg){ inObj->setProjectionLinearTolerance( inArg ); }
PxReal getPxD6Joint_ProjectionLinearTolerance( const PxD6Joint* inObj ) { return inObj->getProjectionLinearTolerance(); }
void setPxD6Joint_ProjectionAngularTolerance( PxD6Joint* inObj, PxReal inArg){ inObj->setProjectionAngularTolerance( inArg ); }
PxReal getPxD6Joint_ProjectionAngularTolerance( const PxD6Joint* inObj ) { return inObj->getProjectionAngularTolerance(); }
const char * getPxD6Joint_ConcreteTypeName( const PxD6Joint* inObj ) { return inObj->getConcreteTypeName(); }
 PxD6JointGeneratedInfo::PxD6JointGeneratedInfo()
	: Motion( "Motion", setPxD6Joint_Motion, getPxD6Joint_Motion)
	, LinearLimit( "LinearLimit", setPxD6Joint_LinearLimit, getPxD6Joint_LinearLimit)
	, TwistLimit( "TwistLimit", setPxD6Joint_TwistLimit, getPxD6Joint_TwistLimit)
	, SwingLimit( "SwingLimit", setPxD6Joint_SwingLimit, getPxD6Joint_SwingLimit)
	, Drive( "Drive", setPxD6Joint_Drive, getPxD6Joint_Drive)
	, DrivePosition( "DrivePosition", setPxD6Joint_DrivePosition, getPxD6Joint_DrivePosition)
	, DriveVelocity( "DriveVelocity", "linear", "angular", setPxD6Joint_DriveVelocity, getPxD6Joint_DriveVelocity)
	, ProjectionLinearTolerance( "ProjectionLinearTolerance", setPxD6Joint_ProjectionLinearTolerance, getPxD6Joint_ProjectionLinearTolerance)
	, ProjectionAngularTolerance( "ProjectionAngularTolerance", setPxD6Joint_ProjectionAngularTolerance, getPxD6Joint_ProjectionAngularTolerance)
	, ConcreteTypeName( "ConcreteTypeName", getPxD6Joint_ConcreteTypeName)
{}
 PxD6JointGeneratedValues::PxD6JointGeneratedValues( const PxD6Joint* inSource )
		:PxJointGeneratedValues( inSource )
		,LinearLimit( getPxD6Joint_LinearLimit( inSource ) )
		,TwistLimit( getPxD6Joint_TwistLimit( inSource ) )
		,SwingLimit( getPxD6Joint_SwingLimit( inSource ) )
		,DrivePosition( getPxD6Joint_DrivePosition( inSource ) )
		,ProjectionLinearTolerance( getPxD6Joint_ProjectionLinearTolerance( inSource ) )
		,ProjectionAngularTolerance( getPxD6Joint_ProjectionAngularTolerance( inSource ) )
		,ConcreteTypeName( getPxD6Joint_ConcreteTypeName( inSource ) )
{
		for ( PxU32 idx = 0; idx < static_cast<PxU32>( physx::PxD6Axis::eCOUNT ); ++idx )
		Motion[idx] = getPxD6Joint_Motion( inSource, static_cast< PxD6Axis::Enum >( idx ) );
		for ( PxU32 idx = 0; idx < static_cast<PxU32>( physx::PxD6Drive::eCOUNT ); ++idx )
		Drive[idx] = getPxD6Joint_Drive( inSource, static_cast< PxD6Drive::Enum >( idx ) );
	getPxD6Joint_DriveVelocity( inSource, DriveVelocity[0], DriveVelocity[1] );
}
inline PxReal getPxJointLimitParametersRestitution( const PxJointLimitParameters* inOwner ) { return inOwner->restitution; }
inline void setPxJointLimitParametersRestitution( PxJointLimitParameters* inOwner, PxReal inData) { inOwner->restitution = inData; }
inline PxReal getPxJointLimitParametersSpring( const PxJointLimitParameters* inOwner ) { return inOwner->spring; }
inline void setPxJointLimitParametersSpring( PxJointLimitParameters* inOwner, PxReal inData) { inOwner->spring = inData; }
inline PxReal getPxJointLimitParametersDamping( const PxJointLimitParameters* inOwner ) { return inOwner->damping; }
inline void setPxJointLimitParametersDamping( PxJointLimitParameters* inOwner, PxReal inData) { inOwner->damping = inData; }
inline PxReal getPxJointLimitParametersContactDistance( const PxJointLimitParameters* inOwner ) { return inOwner->contactDistance; }
inline void setPxJointLimitParametersContactDistance( PxJointLimitParameters* inOwner, PxReal inData) { inOwner->contactDistance = inData; }
 PxJointLimitParametersGeneratedInfo::PxJointLimitParametersGeneratedInfo()
	: Restitution( "Restitution", setPxJointLimitParametersRestitution, getPxJointLimitParametersRestitution )
	, Spring( "Spring", setPxJointLimitParametersSpring, getPxJointLimitParametersSpring )
	, Damping( "Damping", setPxJointLimitParametersDamping, getPxJointLimitParametersDamping )
	, ContactDistance( "ContactDistance", setPxJointLimitParametersContactDistance, getPxJointLimitParametersContactDistance )
{}
 PxJointLimitParametersGeneratedValues::PxJointLimitParametersGeneratedValues( const PxJointLimitParameters* inSource )
		:Restitution( inSource->restitution )
		,Spring( inSource->spring )
		,Damping( inSource->damping )
		,ContactDistance( inSource->contactDistance )
{
}
inline PxReal getPxJointLimitValue( const PxJointLimit* inOwner ) { return inOwner->value; }
inline void setPxJointLimitValue( PxJointLimit* inOwner, PxReal inData) { inOwner->value = inData; }
 PxJointLimitGeneratedInfo::PxJointLimitGeneratedInfo()
	: Value( "Value", setPxJointLimitValue, getPxJointLimitValue )
{}
 PxJointLimitGeneratedValues::PxJointLimitGeneratedValues( const PxJointLimit* inSource )
		:PxJointLimitParametersGeneratedValues( inSource )
		,Value( inSource->value )
{
}
inline PxReal getPxJointLimitPairUpper( const PxJointLimitPair* inOwner ) { return inOwner->upper; }
inline void setPxJointLimitPairUpper( PxJointLimitPair* inOwner, PxReal inData) { inOwner->upper = inData; }
inline PxReal getPxJointLimitPairLower( const PxJointLimitPair* inOwner ) { return inOwner->lower; }
inline void setPxJointLimitPairLower( PxJointLimitPair* inOwner, PxReal inData) { inOwner->lower = inData; }
 PxJointLimitPairGeneratedInfo::PxJointLimitPairGeneratedInfo()
	: Upper( "Upper", setPxJointLimitPairUpper, getPxJointLimitPairUpper )
	, Lower( "Lower", setPxJointLimitPairLower, getPxJointLimitPairLower )
{}
 PxJointLimitPairGeneratedValues::PxJointLimitPairGeneratedValues( const PxJointLimitPair* inSource )
		:PxJointLimitParametersGeneratedValues( inSource )
		,Upper( inSource->upper )
		,Lower( inSource->lower )
{
}
inline PxReal getPxJointLimitConeYAngle( const PxJointLimitCone* inOwner ) { return inOwner->yAngle; }
inline void setPxJointLimitConeYAngle( PxJointLimitCone* inOwner, PxReal inData) { inOwner->yAngle = inData; }
inline PxReal getPxJointLimitConeZAngle( const PxJointLimitCone* inOwner ) { return inOwner->zAngle; }
inline void setPxJointLimitConeZAngle( PxJointLimitCone* inOwner, PxReal inData) { inOwner->zAngle = inData; }
 PxJointLimitConeGeneratedInfo::PxJointLimitConeGeneratedInfo()
	: YAngle( "YAngle", setPxJointLimitConeYAngle, getPxJointLimitConeYAngle )
	, ZAngle( "ZAngle", setPxJointLimitConeZAngle, getPxJointLimitConeZAngle )
{}
 PxJointLimitConeGeneratedValues::PxJointLimitConeGeneratedValues( const PxJointLimitCone* inSource )
		:PxJointLimitParametersGeneratedValues( inSource )
		,YAngle( inSource->yAngle )
		,ZAngle( inSource->zAngle )
{
}
inline PxReal getPxD6JointDriveSpring( const PxD6JointDrive* inOwner ) { return inOwner->spring; }
inline void setPxD6JointDriveSpring( PxD6JointDrive* inOwner, PxReal inData) { inOwner->spring = inData; }
inline PxReal getPxD6JointDriveDamping( const PxD6JointDrive* inOwner ) { return inOwner->damping; }
inline void setPxD6JointDriveDamping( PxD6JointDrive* inOwner, PxReal inData) { inOwner->damping = inData; }
inline PxReal getPxD6JointDriveForceLimit( const PxD6JointDrive* inOwner ) { return inOwner->forceLimit; }
inline void setPxD6JointDriveForceLimit( PxD6JointDrive* inOwner, PxReal inData) { inOwner->forceLimit = inData; }
inline PxD6JointDriveFlags getPxD6JointDriveFlags( const PxD6JointDrive* inOwner ) { return inOwner->flags; }
inline void setPxD6JointDriveFlags( PxD6JointDrive* inOwner, PxD6JointDriveFlags inData) { inOwner->flags = inData; }
 PxD6JointDriveGeneratedInfo::PxD6JointDriveGeneratedInfo()
	: Spring( "Spring", setPxD6JointDriveSpring, getPxD6JointDriveSpring )
	, Damping( "Damping", setPxD6JointDriveDamping, getPxD6JointDriveDamping )
	, ForceLimit( "ForceLimit", setPxD6JointDriveForceLimit, getPxD6JointDriveForceLimit )
	, Flags( "Flags", setPxD6JointDriveFlags, getPxD6JointDriveFlags )
{}
 PxD6JointDriveGeneratedValues::PxD6JointDriveGeneratedValues( const PxD6JointDrive* inSource )
		:Spring( inSource->spring )
		,Damping( inSource->damping )
		,ForceLimit( inSource->forceLimit )
		,Flags( inSource->flags )
{
}
