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

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "wallContactInfo.H"

#include "interAdhesion.H"
#include "wallMatInfo.H"

#include "virtualMeshLevel.H"
#include "wallPlaneInfo.H"

using namespace Foam;
//---------------------------------------------------------------------------//
wallContactInfo::wallContactInfo
(
    ibContactClass& cClass,
    ibContactVars& cVars
)
:
ibContactClass_(cClass),
ibContactVars_(cVars)
{
    bodyId_ = cVars.bodyId_;
    const materialInfo& cMatInfo(cClass.getMatInfo());

    const List<string> cntPatches = wallMatInfo::getWallPatches();
    forAll(cntPatches, patchI)
    {
        const materialInfo& wInfo = wallMatInfo::getWallMatInfo()[cntPatches[patchI]];

        scalar adhPot = 0;
        string adhPotKey;

        if(cMatInfo.getMaterial() < wInfo.getMaterial())
        {
            adhPotKey = cMatInfo.getMaterial();
            adhPotKey += "-";
            adhPotKey += wInfo.getMaterial();
        }
        else
        {
            adhPotKey = wInfo.getMaterial();
            adhPotKey += "-";
            adhPotKey += cMatInfo.getMaterial();
        }
        if(interAdhesion::getInterAdhesion().found(adhPotKey))
        {
            adhPot = interAdhesion::getInterAdhesion()[adhPotKey];
        }

        // compute mean model parameters
        scalar aY = 1/((1 - sqr(cMatInfo.getNu()))/cMatInfo.getY()
            + (1 - sqr(wInfo.getNu()))/wInfo.getY());
        scalar aG = 1/(2*(2 - cMatInfo.getNu())*(1 + cMatInfo.getNu())
            /cMatInfo.getY() + 2*(2 - wInfo.getNu())
            *(1 + wInfo.getNu())/wInfo.getY());
        scalar aGammaN = (cMatInfo.getGamma()*wInfo.getGamma())
            /(cMatInfo.getGamma()+wInfo.getGamma()+SMALL);
        scalar aGammat = (cMatInfo.getGamma()*wInfo.getGamma())
            /(cMatInfo.getGamma()+wInfo.getGamma()+SMALL);
        scalar aMu = (cMatInfo.getMu()+wInfo.getMu())/2;
        scalar maxAdhN = cMatInfo.getAdhN() + wInfo.getAdhN() - 2*adhPot;

        wallMeanPars_.insert(
            cntPatches[patchI],
            physicalProperties(aY, aG,  aGammaN, aGammat, aMu, maxAdhN, 0, 0)
        );
    }

    reduceM_ = 0;
}

wallContactInfo::~wallContactInfo()
{
}
//---------------------------------------------------------------------------//
boundBox wallContactInfo::constructBoundBox
(
    boundBox& bodyBB
)
{
    vector minPoint = bodyBB.min();
    vector maxPoint = bodyBB.max();

    for(int i=0;i<3;i++)
    {    
        minPoint[i] = floor(minPoint[i]/virtualMeshLevel::getCharCellSize())*virtualMeshLevel::getCharCellSize();
        maxPoint[i] = ceil(maxPoint[i]/virtualMeshLevel::getCharCellSize())*virtualMeshLevel::getCharCellSize();
    }
    boundBox BB = boundBox(minPoint,maxPoint);
    return BB;
}
//---------------------------------------------------------------------------//
void wallContactInfo::constructSM()
{
    boundBox bodyBB = ibContactClass_.getGeomModel().getBounds();
    
    boundBox BB = constructBoundBox(bodyBB);

    vector cellNVector = vector(
        ceil((BB.span()[0]/virtualMeshLevel::getCharCellSize())),
        ceil((BB.span()[1]/virtualMeshLevel::getCharCellSize())),
        ceil((BB.span()[2]/virtualMeshLevel::getCharCellSize()))
    );

    // vector bBoxShiftVector = (cellNVector*virtualMeshLevel::getCharCellSize() - bodyBB.span())/2;

    // InfoH << DEM_Info << " -- SM bBoxShiftVector  : " << bBoxShiftVector << endl;
    // InfoH << DEM_Info << " -- SM BB.span()        : " << BB.span() << endl;

    // InfoH << DEM_Info << " -- SM BoundBox     : " << BB.min() <<", " << BB.max() << endl;
    // InfoH << DEM_Info << " -- SM bodyBB       : " << bodyBB.min() <<", " << bodyBB.max() << endl;
    // InfoH << DEM_Info << " -- SM cellNVector  : " << cellNVector << endl;
    // InfoH << DEM_Info << " -- SM charCellSize : " << virtualMeshLevel::getCharCellSize() << endl;
    const scalar charCellSize(virtualMeshLevel::getCharCellSize());

    SM_.reset(new spectatorMesh(cellNVector,BB,charCellSize));
}
//---------------------------------------------------------------------------//
bool wallContactInfo::isInsidePlane(
    vector checkedPoint,
    string wallName
)
{
    const vector& normalVector = wallPlaneInfo::getWallPlaneInfo()[wallName][0];
    const vector& centerPoint = wallPlaneInfo::getWallPlaneInfo()[wallName][1];
    vector testVector = checkedPoint - centerPoint;
    if((testVector & normalVector) < 0)
    {
        return true;
    }
    return false;
}
//---------------------------------------------------------------------------//
bool wallContactInfo::detectWallContact()
{
    contactPatches_.clear();
    clearOldContact();
    boundBox bodyBB = ibContactClass_.getGeomModel().getBounds();
    pointField bBpoints = bodyBB.points();
    const List<string> wallPatches = wallMatInfo::getWallPatches();
    List<bool> isPatchInContact;
    isPatchInContact.setSize(wallPatches.size(),false);
    bool bBContact(false);
    forAll (bBpoints,bBP)
    {
        forAll(wallPatches,wP)
        {
            if(!isInsidePlane(bBpoints[bBP],wallPatches[wP]))
            {
                isPatchInContact[wP] = true;
                bBContact = true;
            }
        }
    }

    forAll(isPatchInContact,wP)
    {
        if(isPatchInContact[wP])
        {
            contactPatches_.append(wallPatches[wP]);
        }
    }

    return bBContact;

}    
//---------------------------------------------------------------------------//
void wallContactInfo::findContactAreas()
{
    autoPtr<DynamicVectorList> contactSTLPoints(
        new DynamicVectorList);

    // InfoH << DEM_Info <<" -- SM bbBoxIsInContact " << endl;        
    // InfoH << DEM_Info <<" -- SM contactPatches_.size() "<< contactPatches_.size() << endl;        
    pointField bodyPoints = ibContactClass_.getGeomModel().getSTLBodyPoints();
    forAll(bodyPoints,bP)
    {
        forAll(contactPatches_,wP)
        {
            if(!isInsidePlane(bodyPoints[bP],contactPatches_[wP]))
            {
                contactSTLPoints().append(bodyPoints[bP]);
            }

        }
    }
    // InfoH << DEM_Info <<" -- SM contactSTLPoints.size() " << contactSTLPoints().size() << endl;

    if(contactSTLPoints().size() > SMALL)
    {
        constructSM();
        // InfoH << DEM_Info << " -- ConstructSM - Is Allive : #1 " << endl;
        List<DynamicList<vector>> possibleSMContact = detectPossibleSMContact(contactSTLPoints(),contactPatches_);
        // InfoH << DEM_Info << " -- ConstructSM - Is Allive : #2 " << endl;
        InfoH << DEM_Info <<" -- SM possibleSCList.size() " << possibleSMContact.size() << endl;
        InfoH << DEM_Info <<" -- SM contactPatches_.size() " << contactPatches_.size() << endl;

        
        forAll(possibleSMContact,SC)
        {
            label emptyCells(0);
            pointField overallContactPoints;
            List<string> contactPatches;
            // InfoH << DEM_Info <<" -- SM possibleSCList["<<SC<<"]"<< possibleSMContact[SC].size() << endl;
            List<Tuple2<point,boundBox>> sMExportList;
            List<Tuple2<point,boundBox>> sMPlaneList;
            List<Tuple2<point,boundBox>> sMInternal;
            
            boundBox cBbox;
            forAll(possibleSMContact[SC],item)
            {
                Tuple2<point,boundBox> sMExport;
                forAll(contactPatches_,wP)
                {
                    boundBox sMEBB = correctSMBBforWall(SM_().getElementBB(possibleSMContact[SC][item]),contactPatches_[wP]);
                    if(sMEBB.minDim() != 0)
                    {
                        sMExport.first() = SM_()[possibleSMContact[SC][item]].initPoint;
                        sMExport.second() = sMEBB;

                        overallContactPoints.append(sMEBB.points());

                        if(!SM_()[possibleSMContact[SC][item]].isInternal)
                        {
                            sMExportList.append(sMExport);
                        }
                        else
                        {
                            sMInternal.append(sMExport);
                        }                    
                        break;
                    }
                    else
                    {
                        emptyCells++;
                    }
                }
            }
            forAll(contactPatches_,wP)
            {

                boundBox sCBBox = getSCBBox(possibleSMContact[SC]);
                List<bool> isPathechInSC;
                isPathechInSC.setSize(contactPatches_.size(),false);
                pointField sCBBPoints = sCBBox.points();
                forAll(sCBBPoints,bBP)
                {
                    if(!isInsidePlane(sCBBPoints[bBP],contactPatches_[wP]))
                    {
                        isPathechInSC[wP]+= true;
                    }
                }
                if(isPathechInSC[wP])
                {
                    cBbox = boundBox(overallContactPoints,false);
                    boundBox planeBox = contactPlaneBBox(constructVMBox(cBbox,contactPatches_[wP]),contactPatches_[wP]);

                    contactPatches.append(contactPatches_[wP]);

                    forAll(contactPatches_,wP2)
                    {
                        if(wP != wP2)
                        {
                            forAll(sCBBPoints,bBP)
                            {
                                if(!isInsidePlane(sCBBPoints[bBP],contactPatches_[wP2]))
                                {
                                    isPathechInSC[wP2]+= true;
                                }
                            }

                            if(isPathechInSC[wP2])
                            {
                                pointField bbPoints = planeBox.points();
                                pointField newBBPoints;
                                forAll(bbPoints,bBP)
                                {
                                    if(!isInsidePlane(bbPoints[bBP],contactPatches_[wP2]))
                                    {
                                        newBBPoints.append(getPlanePoint(bbPoints[bBP],contactPatches_[wP2]));
                                    }
                                }
                                planeBox = boundBox(newBBPoints,false);
                            }

                        }
                    }
                    // InfoH << DEM_Info << "-- SM PlaneBox "<< planeBox << endl;
                    Tuple2<point,boundBox> sMPlane(planeBox.midpoint(),planeBox);
                    sMPlaneList.append(sMPlane);                    
                }

            }                
            // InfoH << DEM_Info <<" -- SM planeBoxList.size() : " << sMPlaneList.size() << endl;
            // InfoH << DEM_Info <<" -- SM sMExportList.size() : " << sMExportList.size() << endl;
            // InfoH << DEM_Info <<" -- SM sMInternal.size()   : " << sMInternal.size() << endl;
            // InfoH << DEM_Info <<" -- SM emptyCells.size()   : " << emptyCells << endl;
            setNewSubContact(sMExportList,sMPlaneList,contactPatches,sMInternal,cBbox);
        }               
    }
}
//---------------------------------------------------------------------------//
List<DynamicList<vector>> wallContactInfo::detectPossibleSMContact
(
    DynamicVectorList& contactPoints,
    List<string>& contactPatches
)
{
    vectorHashSet possibleContactElements;
    vectorHashSet checkedOctreeFaces;
    List<DynamicVectorList> baseSubContactList;
    label iterMax(SM_().matrixSize_[0]*SM_().matrixSize_[1]*SM_().matrixSize_[2]);

    autoPtr<DynamicVectorList> contactElements(
        new DynamicVectorList);

    // InfoH << DEM_Info <<"  -- ConstructSM - Is Allive : #1.1 " << endl;

    forAll(contactPoints,cP)
    {
        vector elementIndex(SM_().getSMCentroidIndex(contactPoints[cP]));
        if(!possibleContactElements.found(elementIndex))
        {
            possibleContactElements.insert(elementIndex);
            SM_()[elementIndex].stlFound = true;
            SM_()[elementIndex].initPointSet = true;
            SM_()[elementIndex].initPoint = contactPoints[cP];
            contactElements().append(elementIndex);
        }
    }
    // InfoH << DEM_Info <<"  -- ConstructSM - Is Allive : #1.2 " << endl;
    autoPtr<DynamicVectorList> nextToCheck(
        new DynamicVectorList);

    autoPtr<DynamicVectorList> auxToCheck(
            new DynamicVectorList);

    autoPtr<DynamicVectorList> subContactElements(
            new DynamicVectorList);

    label iterCount2(0);
    // InfoH << DEM_Info <<"  -- ConstructSM - Is Allive : #1.3 " << endl;
    while((contactElements().size() > SMALL) and iterCount2++ < iterMax)
    {
        label iterCount(0);
        nextToCheck().clear();
        subContactElements().clear();

        nextToCheck().append(contactElements()[0]);
        // nextToCheck().append(edgeFaceNeighbours(mesh,commonFaces[0]));
        while ((nextToCheck().size() > 0) and iterCount++ < iterMax)
        {  	
            auxToCheck().clear();
            forAll(nextToCheck(),cE)
            {
                if (SM_()[nextToCheck()[cE]].toCheck)
                {
                    checkSMElement(nextToCheck()[cE],contactPatches);
                    bool isInMesh(false);
                    bool isInBody(false);
                    bool allInMesh(false);
                    bool allInBody(false);
                    checkElement(nextToCheck()[cE],isInMesh,isInBody,allInMesh,allInBody);
                    if(isInBody && isInMesh)
                    {
                        if(!allInMesh)
                        {   
                            subContactElements().append(nextToCheck()[cE]);
                            auxToCheck().append(SM_().faceNeighbourElements(nextToCheck()[cE]));
                            checkedOctreeFaces.insert(nextToCheck()[cE]);
                            if(allInBody)
                            {
                                SM_()[nextToCheck()[cE]].isInternal = allInBody;
                            }
                        }
                    }                    
                    if((isInBody && !isInMesh)|| SM_()[nextToCheck()[cE]].stlFound)
                    {
                        subContactElements().append(nextToCheck()[cE]);
                        checkedOctreeFaces.insert(nextToCheck()[cE]);
                        auxToCheck().append(SM_().faceNeighbourElements(nextToCheck()[cE]));
                        if(allInBody)
                        {
                            SM_()[nextToCheck()[cE]].isInternal = allInBody;
                        }
                    }
                }
            }
            const autoPtr<DynamicVectorList> helpPtr(nextToCheck.ptr());
            nextToCheck.set(auxToCheck.ptr());
            auxToCheck = helpPtr;
        }
        if(subContactElements().size() > 0)
        {
            baseSubContactList.append(subContactElements());
        }
        for (auto it = contactElements().begin(); it != contactElements().end();)
        {
            if (checkedOctreeFaces.found(*it))
            {
                it = contactElements().erase(it);
            }
            else
            {
                ++it;
            }
        }
    }
    return baseSubContactList;
}
//---------------------------------------------------------------------------//
boundBox wallContactInfo::getSCBBox
(     
    DynamicVectorList& subContactAreas
)
{
    pointField contactArea;
    forAll(subContactAreas,sCE)
    {
        contactArea.append(SM_()[subContactAreas[sCE]].center);
        List<vector> vertexList = SM_().elementVertexIndexies(subContactAreas[sCE]);
        forAll(vertexList,vL)
        {
            contactArea.append(SM_()(vertexList[vL]).center);
        }
    }
    return boundBox(contactArea,false);
}
//---------------------------------------------------------------------------//
boundBox wallContactInfo::constructVMBox
(     
    boundBox& baseContactAreaBB,
    string& wallName
)
{
    // InfoH << DEM_Info << " --SM initialBB Box : " << baseContactAreaBB << endl;
    pointField bbPoints = baseContactAreaBB.points();
    pointField newBBPoints;
    forAll(bbPoints,bP)
    {
        if(isInsidePlane(bbPoints[bP],wallName))
        {
            newBBPoints.append(getPlanePoint(bbPoints[bP],wallName));
        }
        else
        {
           newBBPoints.append(bbPoints[bP]);
        }
    }
    return boundBox(newBBPoints,false);
}
//---------------------------------------------------------------------------//
vector wallContactInfo::getPlanePoint
(
    vector pointInDomain,
    string wallName
)
{
    //contactPlane Data
    const vector& normalVector = wallPlaneInfo::getWallPlaneInfo()[wallName][0];
    const vector& centerPoint = wallPlaneInfo::getWallPlaneInfo()[wallName][1];
    //contactPlane Data
    
    scalar fracI = normalVector & (centerPoint - pointInDomain);
    scalar fracII = magSqr(normalVector);
    // scalar fracI = normalVector[0]*(centerPoint[0]-pointInDomain[0])+
    //                 normalVector[1]*(centerPoint[1]-pointInDomain[1])+
    //                 normalVector[2]*(centerPoint[2]-pointInDomain[2]);
    // scalar fracII = pow(normalVector[0],2) + pow(normalVector[1],2) + pow(normalVector[2],2);

    return pointInDomain + (fracI/fracII)*normalVector;  
}
//---------------------------------------------------------------------------//
void wallContactInfo::checkSMElement
(     
    vector& index,
    List<string> contactPatches
)
{
    label iC(0);
    if(SM_()[index].toCheck)
    {
        point centroidPoint = SM_()[index].center;
        SM_()[index].toCheck = false;
        SM_()[index].isCBody = ibContactClass_.getGeomModel().pointInside(centroidPoint);
        bool isMeshLocC(true);
        while(iC < contactPatches.size() and !SM_()[index].isMesh)
        {
            
            isMeshLocC *= isInsidePlane(centroidPoint,contactPatches[iC]);
            iC++;
        }

        SM_()[index].isMesh = isMeshLocC;

        const List<vector> verticesLabels = SM_().elementVertexIndexies(index);
        forAll(verticesLabels, vL)
        {
            iC = 0;       
            if(SM_()(verticesLabels[vL]).toCheck)
            {
                point vertexPoint = SM_()(verticesLabels[vL]).center;
                SM_()(verticesLabels[vL]).toCheck = false;
                SM_()(verticesLabels[vL]).isCBody = ibContactClass_.getGeomModel().pointInside(vertexPoint);
                bool isMeshLocV(true);
                while(iC < contactPatches.size() and !SM_()(verticesLabels[vL]).isMesh)
                {
                    isMeshLocV *= isInsidePlane(vertexPoint,contactPatches[iC]);
                    iC ++;
                }
                SM_()(verticesLabels[vL]).isMesh = isMeshLocV;
            }
        }
    }
}
//---------------------------------------------------------------------------//
void wallContactInfo::checkElement
(
    vector& index,
    bool& inMesh,
    bool& inBody,
    bool& allInMesh,
    bool& allInBody  
)
{
    bool inMeshLocal(false);
    bool inBodyLocal(false);
    bool allInMeshLocal(true);
    bool allInBodyLocal(true);

    inMeshLocal =  inMeshLocal || SM_()[index].isMesh;
    allInMeshLocal *= SM_()[index].isMesh;
    inBodyLocal =  inBodyLocal || SM_()[index].isCBody;
    allInBodyLocal *= SM_()[index].isCBody;
    if(SM_()[index].isCBody && !SM_()[index].initPointSet)
    {
        SM_()[index].initPoint = SM_()[index].center;
        SM_()[index].initPointSet = true;
    }
    const List<vector> verticesLabels = SM_().elementVertexIndexies(index);

    forAll(verticesLabels, vL)
    {     
        inMeshLocal =  inMeshLocal || SM_()(verticesLabels[vL]).isMesh;
        allInMeshLocal *= SM_()(verticesLabels[vL]).isMesh;
        inBodyLocal =  inBodyLocal || SM_()(verticesLabels[vL]).isCBody;
        allInBodyLocal *= SM_()(verticesLabels[vL]).isCBody;
        if(SM_()(verticesLabels[vL]).isCBody && !SM_()[index].isCBody && !SM_()[index].initPointSet)
        {
            SM_()[index].initPoint = SM_()(verticesLabels[vL]).center;
            SM_()[index].initPointSet = true;
        }
    }  

    inMesh = inMeshLocal;
    inBody = inBodyLocal;
    allInMesh = allInMeshLocal;
    allInBody = allInBodyLocal; 
}
//---------------------------------------------------------------------------//
boundBox wallContactInfo::contactPlaneBBox
(
    boundBox contactBoundBox,
    string wallName
)
{
    pointField planeBB;
    const pointField contactBB = contactBoundBox.points();
    forAll(contactBB,bBP)
    {
        planeBB.append(getPlanePoint(contactBB[bBP],wallName));
    }
    return boundBox(planeBB,false);
}
//---------------------------------------------------------------------------//
boundBox wallContactInfo::correctSMBBforWall
(
    boundBox bB,
    string wallName
)
{
    pointField newBB;
    const pointField elementBBPoints = bB.points();
    forAll(elementBBPoints,bBP)
    {
        if(isInsidePlane(elementBBPoints[bBP],wallName))
        {
            newBB.append(getPlanePoint(elementBBPoints[bBP],wallName));
        }
        else
        {
           newBB.append(elementBBPoints[bBP]); 
        }
    }
    return boundBox(newBB,false);
}
//---------------------------------------------------------------------------//
void wallContactInfo::clearOldContact()
{
    subCList_.clear();
}
//---------------------------------------------------------------------------//
void wallContactInfo::setNewSubContact(
    const List<Tuple2<point,boundBox>>& contactBBData,
    const List<Tuple2<point,boundBox>>& planeBBData,
    const List<string>& contactPatches,
    const List<Tuple2<point,boundBox>>& isInternal,
    boundBox BB
)
{
    subCList_.emplace_back(std::make_shared<wallSubContactInfo>
    (contactBBData, planeBBData,contactPatches,isInternal,wallMeanPars_,BB,bodyId_));

}
//---------------------------------------------------------------------------//
void wallContactInfo::registerSubContactList(DynamicList<wallSubContactInfo*>& wallContactList)
{
    for(auto sC : subCList_)
    {
        wallContactList.append(sC.get());
    }
}
//---------------------------------------------------------------------------//
