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
    addModelRepeatRandomPosition.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef addModelRepeatRandomPosition_H
#define addModelRepeatRandomPosition_H

#include "addModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class addModelRepeatRandomPosition Declaration
\*---------------------------------------------------------------------------*/

class addModelRepeatRandomPosition
:
    public addModel
{

    private:
    typedef DynamicList<label>  DynamicLabelList;
    const dictionary&           addModelDict_;
    word                        addMode_;
    bool                        bodyAdded_;

    const dictionary&   coeffsDict_;

	word			    addDomain_;
    word			    scalingMode_;
    word			    rotationMode_;
	word			    addModeI_;

	const dictionary& addDomainCoeffs_;
	const dictionary& scalingModeCoeffs_;
	const dictionary& rotationModeCoeffs_;
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

	bool 		scaleParticles_;
	scalar 		minScale_;
	scalar		maxScale_;
	scalar 		minScaleFit_;
	scalar		scaleStep_;
	scalar		nTriesBeforeScaling_;

	bool 		rotateParticles_;
	bool		randomAxis_;
    vector 		axisOfRot_;

	int 		bodyAdditionAttemptCounter_;
	int 		scaleCorrectionCounter_;

	bool 		scaleApplication_;
	bool 		scaleRandomApplication_;
	bool 		rescaleRequirement_;
	bool 		succesfulladition_;
	scalar		scalingFactor_;
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

    public:

    // Constructors ------------------------------------------------
    addModelRepeatRandomPosition
    (
        const dictionary& addModelDict,
        const Foam::fvMesh& mesh,
        std::unique_ptr<geomModel> bodyGeomModel,
        List<labelList>& cellPoints
    );

    // Destructors -------------------------------------------------
    ~addModelRepeatRandomPosition();

    // Public methods ----------------------------------------------
    void init();

    bool shouldAddBody(const volScalarField& body);

    std::shared_ptr<geomModel> addBody(
        const volScalarField& body,
        PtrList<immersedBody>& immersedBodies
    );

	scalar returnRandomAngle();

	scalar returnRandomScale();

    vector returnRandomRotationAxis();

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

    bool checkIfAdd()
    {
        return true;
    }

    word getAddMode()
    {
        return addMode_;
    }

    bool getBodyAdded()
    {
        return bodyAdded_;
    }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
