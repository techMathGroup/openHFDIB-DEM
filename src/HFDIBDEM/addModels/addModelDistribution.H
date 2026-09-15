/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  | |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___/  |_____|_|  |_|
      | |                     H ybrid F ictitious D omain - I mmersed B oundary
      |_|                                        and D iscrete E lement M ethod
-------------------------------------------------------------------------------
License

    openHFDIB-DEM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (Version 3) as published
    by the Free Software Foundation.

    openHFDIB-DEM is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with openHFDIB-DEM. If not, see <http://www.gnu.org/licenses/>.

InNamespace
    Foam

Description
    class for addition and removal of immersedBodies

SourceFiles
    addModelDistribution.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef addModelDistribution_H
#define addModelDistribution_H

#include "addModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class addModelDistribution Declaration
\*---------------------------------------------------------------------------*/

class addModelDistribution
:
    public addModel
{

    private:
    typedef DynamicList<label>  DynamicLabelList;
    typedef DynamicList<scalar>  DynamicScalarList;
    const dictionary&           addModelDict_;
    word                        addMode_;
    bool                        bodyAdded_;

    IOdictionary           distributionDict_;
    const scalarList       distribution_;
    const scalarList       particleSize_;
    const scalar       convertToMeters_;

    DynamicScalarList  addedParticlesSize_;
    scalar               volumeOfAddedBodies_;

    const dictionary&   coeffsDict_;
    const scalar       stlBaseSize_;

	word			    addDomain_;
	word			    addModeI_;

	const dictionary& addDomainCoeffs_;
	const dictionary& addModeICoeffs_;

    label       useNTimes_;
    scalar      timeBetweenUsage_;
    label		partPerAdd_;
	scalar		fieldValue_;
    label		addedOnTimeLevel_;
    label		partPerAddTemp_;

	word		zoneName_;
    vector		minBound_;
    vector		maxBound_;

	int 		bodyAdditionAttemptCounter_;
	int 		scaleCorrectionCounter_;

	bool 		succesfulladition_;
	bool		restartPartCountTemp_;
	bool		reapeatedAddition_;
	bool		firstTimeRunning_;
	bool 		cellZoneActive_;
	bool		boundBoxActive_;
	Field<label>octreeField_;
	bool		timeBased_;
	bool		fieldBased_;
	scalar		fieldCurrentValue_;

    List<DynamicLabelList>  cellsInBoundBox_;
    List<pointField>        cellZonePoints_;
    boundBox                cellZoneBounds_;

    bool        allActiveCellsInMesh_;

    Random          randGen_;

    protected:


    public:

    // Constructors ------------------------------------------------
    addModelDistribution
    (
        const dictionary& addModelDict,
        const Foam::fvMesh& mesh,
        std::unique_ptr<geomModel> bodyGeomModel,
        List<labelList>& cellPoints
    );

    // Destructors -------------------------------------------------
    ~addModelDistribution();

    // Public methods ----------------------------------------------
    void init();

    bool shouldAddBody(const volScalarField& body);

    std::shared_ptr<geomModel> addBody(
        const volScalarField& body,
        PtrList<immersedBody>& immersedBodies
        );

    Tuple2<label, scalar> returnScaleFactor();

	void initializeCellZone();

	void updateCellZoneBoundBox();

    void initializeBoundBox();

    void recreateBoundBox();

	labelList getBBoxCellsByOctTree
    (
        label cellToCheck,
        bool& insideBB,
        vector& bBoxMin,
        vector& bBoxMax,
        List<DynamicLabelList>& bBoxCells
    );

	scalar checkLambdaFraction(const volScalarField& body);

    bool checkIfAdd() {return true;}

    word getAddMode()
    {
        return addMode_;
    }

    bool getBodyAdded()
    {
        return bodyAdded_;
    }

    scalar returnRandomAngle();

    vector returnRandomRotationAxis();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
