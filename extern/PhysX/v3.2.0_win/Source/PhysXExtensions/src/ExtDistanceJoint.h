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


#ifndef NP_DISTANCEJOINTCONSTRAINT_H
#define NP_DISTANCEJOINTCONSTRAINT_H

#include "PsUserAllocated.h"
#include "ExtJoint.h"
#include "PxDistanceJoint.h"
#include "PxTolerancesScale.h"

namespace physx
{
namespace Ext
{

	struct DistanceJointData : public JointData
	{
							PxReal					minDistance;
							PxReal					maxDistance;
							PxReal					tolerance;
							PxReal					spring;
							PxReal					damping;

							PxDistanceJointFlags	jointFlags;
		EXPLICIT_PADDING(	PxU16					paddingFromFlags);
	};

typedef Joint<PxDistanceJoint, PxJointType::eDISTANCE> DistanceJointT;

	class DistanceJoint : public DistanceJointT
	{
	public:
// PX_SERIALIZATION
									DistanceJoint(PxRefResolver& v)	: DistanceJointT(v)	{}
									DECLARE_SERIAL_CLASS(DistanceJoint, DistanceJointT)
		virtual		void			exportExtraData(PxSerialStream& stream);
		virtual		char*			importExtraData(char* address, PxU32& totalPadding);
		virtual		bool			resolvePointers(PxRefResolver&, void*);
		static		void			getMetaData(PxSerialStream& stream);
//~PX_SERIALIZATION
		virtual ~DistanceJoint()
		{
			if(getSerialFlags()&PxSerialFlag::eOWNS_MEMORY)
				PX_FREE(mData);
		}

		DistanceJoint(const PxTolerancesScale& scale,
					  PxRigidActor* actor0, const PxTransform& localFrame0, 
					  PxRigidActor* actor1, const PxTransform& localFrame1)
		{
			setSerialType(PxConcreteType::eUSER_DISTANCE_JOINT);
			DistanceJointData* data = reinterpret_cast<DistanceJointData*>(PX_ALLOC(sizeof(DistanceJointData), PX_DEBUG_EXP("DistanceJointData")));
			mData = data;

			initCommonData(*data,actor0, localFrame0, actor1, localFrame1);

			data->spring = 0;
			data->damping = 0;
			data->minDistance = 0;
			data->maxDistance = 0;
			data->tolerance = 0.025f * scale.length;
			data->jointFlags = PxDistanceJointFlag::eMAX_DISTANCE_ENABLED;
		}


		void					setMinDistance(PxReal distance);
		PxReal					getMinDistance()					const;

		void					setMaxDistance(PxReal distance);
		PxReal					getMaxDistance()					const;

		void					setTolerance(PxReal tolerance);
		PxReal					getTolerance()						const;

		void					setSpring(PxReal spring);
		PxReal					getSpring()							const;

		void					setDamping(PxReal damping);
		PxReal					getDamping()						const;
		
		PxDistanceJointFlags	getDistanceJointFlags(void)			const;
		void					setDistanceJointFlags(PxDistanceJointFlags flags);
		void					setDistanceJointFlag(PxDistanceJointFlag::Enum flag, bool value);

		bool					attach(PxPhysics &physics, PxRigidActor* actor0, PxRigidActor* actor1);
	private:

		static PxConstraintShaderTable sShaders;

		PX_FORCE_INLINE DistanceJointData& data() const				
		{	
			return *static_cast<DistanceJointData*>(mData);
		}
	};

} // namespace Ext

namespace Ext
{
	extern "C" PxU32 DistanceJointSolverPrep(Px1DConstraint* constraints,
		PxVec3& body0WorldOffset,
		PxU32 maxConstraints,
		const void* constantBlock,
		const PxTransform& bA2w,
		const PxTransform& bB2w);
}

}

#endif
