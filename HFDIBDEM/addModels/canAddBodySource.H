    // Note: this function is actually common to all the derived classes
    //       but due to the constructor fuckery, it cannot be really 
    //       implemented in the baseClass (AFAIK, please prove me wrong)
    
    // TODO: get this to work with Martas's OCTREE algorithm
    //       (see getBBoxCellsByOctree)
    
    //~ boundBox bounds(bodySurfMesh_.bounds());
    //~ minBoundPoint = bounds.min();
    //~ maxBoundPoint = bounds.max();
    
    //~ scalar inflFact(2*sqrt(mesh_.magSf()[0]));
    
    //~ vector expMinBBox = minBoundPoint_ - vector::one*inflFact;
    //~ vector expMaxBBox = maxBoundPoint_ + vector::one*inflFact;
    
    //~ octreeField_ *= 0;
    //~ bool insideBB(false);
    //~ List<DynamicLabelList> bBoxCells(Pstream::nProcs());
    
    const triSurface ibTemp(bodySurfMesh);
    triSurfaceSearch ibTriSurfSearch( ibTemp );
    const pointField& pp = mesh_.points();
    // Return true only if all cells that will be cover by new body has lambda equal to zero
    forAll (body,cellI)
    {
        //Check if partially or completely inside
        const labelList& vertexLabels = mesh_.cellPoints()[cellI];
        const pointField vertexPoints(pp,vertexLabels);
        boolList vertexesInside = ibTriSurfSearch.calcInside( vertexPoints );
        
        forAll (vertexesInside, verIn)
        {
            if (vertexesInside[verIn]==true && body[cellI] > SMALL)
            {
                return false;
            }
        }
    }
    return true;
