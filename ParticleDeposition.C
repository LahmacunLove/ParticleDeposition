/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ParticleDeposition.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::ParticleDeposition<CloudType>::applyToPatch
(
    const label globalPatchi
) const
{
    forAll(patchIDs_, i)
    {
        if (patchIDs_[i] == globalPatchi)
        {
            return i;
        }
    }

    return -1;
}


template<class CloudType>
void Foam::ParticleDeposition<CloudType>::write()
{
    if (PDRPtr_.valid())
    {
        PDRPtr_->write();
    }
    else
    {
        FatalErrorInFunction
            << "PDRPtr not valid" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleDeposition<CloudType>::ParticleDeposition
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    PDRPtr_(nullptr),
    patchIDs_()
{
    const wordList allPatchNames(owner.mesh().boundaryMesh().names());
    const wordReList patchNames(this->coeffDict().lookup("patches"));

    labelHashSet uniqIds;
    for (const wordRe& re : patchNames)
    {
        labelList ids = findStrings(re, allPatchNames);

        if (ids.empty())
        {
            WarningInFunction
                << "Cannot find any patch names matching " << re
                << endl;
        }

        uniqIds.insert(ids);
    }

    patchIDs_ = uniqIds.sortedToc();

    // trigger creation of the PDR field
    preEvolve();
}


template<class CloudType>
Foam::ParticleDeposition<CloudType>::ParticleDeposition
(
    const ParticleDeposition<CloudType>& pe
)
:
    CloudFunctionObject<CloudType>(pe),
    PDRPtr_(nullptr),
    patchIDs_(pe.patchIDs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleDeposition<CloudType>::preEvolve()
{
    if (PDRPtr_.valid())
    {
        PDRPtr_->primitiveFieldRef() = 0.0;
    }
    else
    {
        const fvMesh& mesh = this->owner().mesh();

        PDRPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "PDR",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimVolume, Zero)
            )
        );
    }
}


template<class CloudType>
void Foam::ParticleDeposition<CloudType>::postPatch
(
    const parcelType& p,
    const polyPatch& pp,
    bool&
)
{
    const label patchi = pp.index();

    const label localPatchi = applyToPatch(patchi);

    if (localPatchi != -1)
    {
        vector nw;
        vector Up;

        // patch-normal direction
        this->owner().patchData(p, pp, nw, Up);

        // particle velocity relative to patch
        const vector& U = p.U() - Up;

        // quick reject if particle travelling away from the patch
        if ((nw & U) < 0)
        {
            return;
        }

        const scalar coeff = p.nParticle()

        const label patchFacei = pp.whichFace(p.face());
        scalar& PDR = PDRPtr_->boundaryFieldRef()[patchi][patchFacei];

        PDR += coeff;
    }
}


// ************************************************************************* //
