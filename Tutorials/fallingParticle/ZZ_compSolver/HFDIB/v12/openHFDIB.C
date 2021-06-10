/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________
                       | | | ||  ___|  _  \_   _| ___ \     H ybrid
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /     F ictitious
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \     D omain
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /     I mmersed
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/      B oundary
      | |
      |_|
-------------------------------------------------------------------------------
License

    openHFDIB is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

Contributors
    Federico Municchi (2016)
\*---------------------------------------------------------------------------*/
#include "openHFDIB.H"
#include "polyMesh.H"
#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"

#include "interpolationCellPoint.H"
#include "interpolationCell.H"

#define ORDER 2

using namespace Foam;

// Q (MI)   : why are there doubled includes?
// Note (MI): removed the doubled includes and STD includes

//---------------------------------------------------------------------------//
openHFDIB::openHFDIB(const Foam::dynamicFvMesh& mesh )
:
mesh_(mesh),
HFDIBDict_
(
    IOobject
    (
        "HFDIBDict",
        "constant",
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
transportProperties_
(
    IOobject
    (
        "transportProperties",
        "constant",
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
)
{
}
//---------------------------------------------------------------------------//
openHFDIB::~openHFDIB()
{
    //~ forAll(immersedBodies_,bodyId)
    //~ {
        //~ delete immersedBodies_[bodyId];
        //~ immersedBodies_[bodyId] = NULL;
    //~ }
    
    immersedBodies_.clear();
    //~ prtContactLst_.clear();
    //~ wallContactLst_.clear();
    //~ pairingLst_.clear();
}
//---------------------------------------------------------------------------//
void openHFDIB::initialize()
{

    wordList stlNames( HFDIBDict_.lookup("stlNames") );
    HFDIBinterpDict_ = HFDIBDict_.subDict("interpolationSchemes");
    //Generate immersed objects
    immersedBodies_.setSize(stlNames.size());
    forAll(stlNames,nameI)
    {
        word stlName(stlNames[nameI]);
        Info << "Loading stl file: " << stlName << endl;
        
        immersedBodies_.set
        (
            nameI,
            new immersedBody
            (
                stlName,
                mesh_,
                HFDIBDict_,
                transportProperties_,
                nameI
            )
        );
    }
}
//---------------------------------------------------------------------------//
void openHFDIB::update
(
    volScalarField& body,
    volVectorField& f
)
{    
    //Contact variables (local to this subroutine)
    // - particles undergoing contact with other particles
    DynamicLabelList    prtContactLst;
    // - particles undergoing contact with walls
    DynamicLabelList    wallContactLst;
    // - which particle interacts with which
    DynamicList<DynamicLabelList> pairingLst;
    

    forAll(immersedBodies_,bodyId)
    {
        // move the body
        immersedBodies_[bodyId].updateBodyField(body,f);
        
        // check for contact
        if (immersedBodies_[bodyId].checkPrtContact())
        {
            prtContactLst.append(immersedBodies_[bodyId].getBodyId());
        }
        else
        {
            Info << "doing nothing" << endl;
        }
        if (immersedBodies_[bodyId].checkWallContact()) 
        {
            wallContactLst.append(immersedBodies_[bodyId].getBodyId());
        }
    }
    
    Info << "Bodies undergoing inter-particle contact:" << endl;
    Info << prtContactLst << endl;
    
    // compute the contact forces acting on each body
    forAll(prtContactLst,lstI)
    {// outer loop through the particles
        // Note (MI): in this loop, I construct the paringList
        
        // current body indentification
        label cBodyId(prtContactLst[lstI]);
        
        // reference to current body
        immersedBody& cBody(immersedBodies_[cBodyId]);
        // referemce to the cells of the current body that are in contact
        // with other particles
        const DynamicLabelList& cPrtContactCells(cBody.getPrtContactCellList());
        
        // get data from the current body necessary for computation of
        // contact forces
        vector cCoM(cBody.getCoM());        //center of mass
        vector cVel(cBody.getVel());        //linear velocity
        scalar cKN(cBody.getKN());          //elastic normal stiffness
        scalar cGammaN(cBody.getGammaN());  //normal viscosity
        
        DynamicLabelList cBodyContacts;
        
        for (label lstII=lstI+1;lstII<prtContactLst.size();lstII++)
        {//inner loop, check with which particles has the current one common points
            label tBodyId(prtContactLst[lstII]);
            
            immersedBody& tBody(immersedBodies_[tBodyId]);
            const DynamicLabelList& tPrtContactCells(tBody.getPrtContactCellList());
            
            // identify the common cells between the current and tested body
            DynamicLabelList commonCells;
            forAll(cPrtContactCells,cContCellI)
            {
                forAll(tPrtContactCells,tContCellI)
                {
                    if (cPrtContactCells[cContCellI] == tPrtContactCells[tContCellI]) commonCells.append(cPrtContactCells[cContCellI]);
                }
            }
            
            // Note (MI): at this moment, I should have available cells
            //            in contact for both the bodies
            //         => I should be able to compute the contact forces
            //            between the two bodies
            if (commonCells.size() > 0)
            {
                cBodyContacts.append(tBodyId);
                Info << "Cells common between body " << cBodyId << " and " << tBodyId << ":" << endl;
                Info << commonCells << endl;
                
                // Note: from now on, I work only with 2 particles
                
                // - get local particle radii
                DynamicScalarList cALst(cBody.getLocPartRad(commonCells));
                DynamicScalarList tALst(tBody.getLocPartRad(commonCells));
                
                // estimate delta
                scalar delta(0.0);
                forAll(commonCells,cCellI)
                {
                    delta+=mesh_.V()[commonCells[cCellI]];
                }
                delta/=commonCells.size();
                delta =Foam::pow(delta,0.33333);
                
                Info << "-- Approximate penetration depth: " << delta << " m" << endl;
                
                // get data from the tested body necessary for computation of
                // contact forces
                vector tCoM(tBody.getCoM());        //center of mass
                vector tVel(tBody.getVel());        //linear velocity
                scalar tKN(tBody.getKN());          //elastic normal stiffness
                scalar tGammaN(tBody.getGammaN());  //normal viscosity
                
                // compute mean model parameters
                scalar aKN(0.5*(cKN+tKN));
                scalar aGammaN(0.5*(cGammaN+tGammaN));
                
                Info << "-- Got the necessary data and updated model parameters" << endl;
                
                //~ scalar redM((cBody.getM()*tBody.getM())/(cBody.getM()+tBody.getM()));
                //~ scalar tC(3.1415/(Foam::sqrt(aKN/(redM+SMALL) - Foam::sqr(0.5*aGammaN/(redM+SMALL)))));
                
                //~ Info << "-- Typical contact response time: " << tC << endl;
                
                // compute the normal force
                vector nVec((cCoM - tCoM)/mag(cCoM - tCoM));
                scalar Vn(-(cVel - tVel) & nVec);
                vector FN((aKN*delta + aGammaN*Vn)*nVec);
                
                Info << "-- Normal force: " << FN << endl;
                
                // add the computed force to the affected bodies
                cBody.updateFAndT(FN,vector::zero);
                tBody.updateFAndT(-FN,vector::zero);
                
                Info << "-- updated F and T in the bodies in contact" << endl;
                
                
                
                //~ Info << "-- Acting normal force" << endl;
                //~ Info << FN << endl;
            }
        }
        
        // update the pairing list - with which bodies is the current body
        // in contact?
        pairingLst.append(cBodyContacts);
    }
    Info << pairingLst << endl;
    
    // move the bodies due to the contact
    forAll(prtContactLst,lstI)
    {
        label cBodyId(prtContactLst[lstI]);
        immersedBody& cBody(immersedBodies_[cBodyId]);
        cBody.postContactBodyUpdate(body);
    }
}
//---------------------------------------------------------------------------//
void openHFDIB::forceReinit
(
    volScalarField& body,
    volVectorField& f
)
{
    forAll(immersedBodies_,bodyId)
    {
        Info << "reseting body " << bodyId << endl;
        immersedBodies_[bodyId].resetIB(body);
        immersedBodies_[bodyId].updateBodyField(body,f);
        Info << "done" << endl;
    }
}
//---------------------------------------------------------------------------//
void openHFDIB::interpolateIB( volVectorField & V
                              ,volVectorField & Vs
                              ,volScalarField & body)
{

    //Create interpolator
    Info << "creating interpolator" << endl;
    autoPtr<interpolation<vector>> interpV =
                   interpolation<vector>::New(HFDIBinterpDict_, V);
    Info << "created interpolator" << endl;
    
    //Reset imposed field
    Vs *= scalar(0);
    
    //Loop over all the immersed bodies
    forAll(immersedBodies_,bodyId)
    {
        //Update imposed field according to body
        immersedBodies_[bodyId].updateVectorField(Vs, V.name(),body);
        
        // Get surface cells, interpolation points and internal cells
        const DynamicLabelList&  surCells  = immersedBodies_[bodyId].getSurfaceCellList();
        const List<DynamicPointList>& intPoints = immersedBodies_[bodyId].getInterpolationPoints();;
        const List<DynamicLabelList>& intCells  = immersedBodies_[bodyId].getInterpolationCells();
        
        //loop over all surface cells
        forAll(surCells,scell)
        {
            label cellI = surCells[scell];
            //Check max order of accuracy
            List<bool> allowedOrder;
            allowedOrder.setSize(ORDER);
            // Note (MI): Ok, i will create a list of bools of size ORDER
            
            for (int intPoint=0;intPoint<ORDER;intPoint++)
            {
                if (intCells[scell][intPoint] == -1)
                {
                    allowedOrder[intPoint] = false;
                }
                else
                {
                    allowedOrder[intPoint] = true;
                }
            }
            // Note (MI): I expanded this for-if for better readability
            
            bool  firstOrder = false;
            bool secondOrder = true;
            bool   zeroOrder = false;
            
            //Check if the second order is possible
            if( allowedOrder[1] == false)
            {
                secondOrder = false;
                firstOrder  = true;
            }
            
            //Check if first order is possible
            if( allowedOrder[0] == false)
            {
                secondOrder = false;
                firstOrder  = false;
                zeroOrder   = true;
            }
            
            //Go for interpolation!
            if(secondOrder)
            {
                //~ Info << "second order" << endl;
            
                vector VP1 =  interpV->interpolate(  intPoints[scell][1],
                                                  intCells[scell][0]
                                                ) - Vs[cellI];
                
                vector VP2 =  interpV->interpolate(  intPoints[scell][2],
                                                  intCells[scell][1]
                                                ) - Vs[cellI];
                
                
                //distance between interpolation points
                scalar res_ = mag(intPoints[scell][2] -intPoints[scell][1]);
                
                //cell center to surface distance
                scalar ds   = res_*(0.5-body[cellI ]) ;

                vector quadCoeff = 1.0/(res_*res_+SMALL) * ( VP2/2.0 - VP1 );
                vector linCoeff  = 1.0/(2.0*res_+SMALL) * ( 4.0*VP1 - VP2 );
                // Note (MI): added +SMALL for the robustness, effects
                //            on solution not tested
                
                Vs[cellI] = quadCoeff*ds*ds + linCoeff*ds + Vs[cellI]  ;
            }
            else if(firstOrder)
            {
                //~ Info << "First order" << endl;
                vector VP1 =  interpV->interpolate(  intPoints[scell][1],
                                                      intCells[scell][0]
                                                    ) - Vs[cellI];
                
                
                //distance between interpolation points
                scalar res_ = mag(intPoints[scell][1] -intPoints[scell][0]);
                
                //cell center to surface distance
                scalar ds   = res_*(0.5-body[cellI]) ;
                
                vector linCoeff = VP1/(res_+SMALL);
                // Note (MI): added +SMALL for robustness, effect on
                //            solution not tested
                
                Vs[cellI] = linCoeff*ds + Vs[cellI];
                //~ Vs[cellI] = -linCoeff*ds - Vs[cellI];
            }
            else if(zeroOrder)
            {
                //Zero order
                Vs[cellI] = body[cellI]*Vs[cellI]  + (1.0-body[cellI])*V[cellI];
            }
            // Note (MI): I reset Vs field prior to this computation,
            //            thus, isn't Vs[cellI] always 0 before assignement?
            //            (from the tests, it would seem so)
            // Note (MI): actually, no. updateVectorField recreates the
            //            Vs field
        
        }
    
    }
}
