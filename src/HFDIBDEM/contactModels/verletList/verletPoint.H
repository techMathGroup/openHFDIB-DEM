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
    class for Verlet point

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef verletPoint_H
#define verletPoint_H

#include "dictionary.H"
#include "IDLList.H"
#include "fvCFD.H"
#include "geomModel.H"

#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class verletPoint Declaration
\*---------------------------------------------------------------------------*/

class verletBox;

class verletPoint
{
    private:

    label bodyID_;
    std::shared_ptr<verletBox> parentBox_;

    point& vPoint_;
    const bool isMin_;
    bool isStatic_;

    public:

    // Constructors ------------------------------------------------
    verletPoint
    (
        label bodyID,
        std::shared_ptr<verletBox> parentBox,
        point& vPoint,
        const bool isMin,
        bool isStatic
    )
    :
    bodyID_(bodyID),
    parentBox_(parentBox),
    vPoint_(vPoint),
    isMin_(isMin),
    isStatic_(isStatic)
    {};

    // Destructors -------------------------------------------------
    virtual ~verletPoint(){};

    bool isMin()
    {
        return isMin_;
    }

    point& getPoint()
    {
        return vPoint_;
    }

    label getBodyId()
    {
        return bodyID_;
    }

    bool getIsStatic()
    {
        return isStatic_;
    }

    std::shared_ptr<verletBox> getParentBox()
    {
        return parentBox_;
    }
};

class verletBox
: public std::enable_shared_from_this<verletBox>
{
    private:

    label bodyID_;
    std::shared_ptr<boundBox> bBox_;
    bool isStatic_;

    std::shared_ptr<verletPoint> minPoint_;
    std::shared_ptr<verletPoint> maxPoint_;

    // Constructors ------------------------------------------------
    verletBox
    (
        label bodyID,
        std::shared_ptr<boundBox> bBox,
        bool isStatic
    )
    :
    bodyID_(bodyID),
    bBox_(bBox),
    isStatic_(isStatic)
    {};

    public:

    std::shared_ptr<verletBox> getptr()
    {
        return shared_from_this();
    }
    // No public constructor, only a factory function,
    // so there's no way to have getptr return nullptr.
    [[nodiscard]] static std::shared_ptr<verletBox> create
    (
        label bodyID,
        std::shared_ptr<boundBox> bBox,
        bool isStatic
    )
    {
        // Not using std::make_shared<Best> because the c'tor is private.
        return std::shared_ptr<verletBox>(new verletBox(bodyID, bBox, isStatic));
    }

    std::shared_ptr<boundBox>& getBBox()
    {
        return bBox_;
    }

    label getBodyId()
    {
        return bodyID_;
    }

    void setVerletPoints()
    {
        minPoint_ = std::make_shared<verletPoint>(bodyID_, shared_from_this(), bBox_->min(), true,isStatic_);
        maxPoint_ = std::make_shared<verletPoint>(bodyID_, shared_from_this(), bBox_->max(), false,isStatic_);
    }

    std::shared_ptr<verletPoint> getMinPoint()
    {
        return minPoint_;
    }

    std::shared_ptr<verletPoint> getMaxPoint()
    {
        return maxPoint_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
