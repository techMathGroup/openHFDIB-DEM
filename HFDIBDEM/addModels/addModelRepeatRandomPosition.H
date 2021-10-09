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

    openHFDIB-DEM is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

Description
    class for addition and removal of immersedBodies

SourceFiles
    addModelRepeatRandomPosition.C

Contributors
    Martin Isoz (2019-*), Martin Šourek (2019-*),
    Ondřej Studeník (2020-*)
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

    autoPtr<geomModel> geomModel_;

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

    label       nGeometricD_;
    Vector<label>      geometricD_;

    Random          randGen_;

    public:

    // Constructors ------------------------------------------------
    addModelRepeatRandomPosition
    (
        const dictionary& addModelDict,
        const Foam::dynamicFvMesh& mesh,
        const Vector<label> geomDir,
        geomModel* bodyGeomModel
    );

    // Destructors -------------------------------------------------
    ~addModelRepeatRandomPosition();

    // Public methods ----------------------------------------------
    void init();

    bool shouldAddBody(const volScalarField& body);

    geomModel* addBody(const volScalarField& body);

	scalar returnRandomAngle();

	scalar returnRandomScale();

    vector returnRandomRotationAxis();

	vector computeCoM(const triSurfaceMesh& bodySurfMesh);

	bool isBodyInBoundBox
    (
        vector minBound_,
        vector maxBound_,
        const triSurfaceMesh& bodySurfMesh
    );

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
