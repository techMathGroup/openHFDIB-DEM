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
    contactModel

Description
    methods to solve particle-particle contact

SourceFiles
    prtContact.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef outputHFDIBDEM_H
#define outputHFDIBDEM_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace outputHFDIBDEM Declaration
\*---------------------------------------------------------------------------*/

namespace outputHFDIBDEM
{
    enum outputType
    {
        basic_Info,
        iB_Info,
        DEM_Info,
        addModel_Info,
        parallelDEM_Info
    };

    class
    InfoHFDIBDEM
    {
        private:
        bool printBasic;
        bool printIb;
        bool printDEM;
        bool printAddModel;
        bool printParallelDEM;
        bool print;
        bool printParDEM;

        HashTable<bool, outputType,Hash<outputType>> outputTable;

        public:
        // Constructors ----------------------------------------------------
        InfoHFDIBDEM()
        {
            printBasic = true;
            printIb = true;
            printDEM = true;
            printAddModel = true;
            printParallelDEM = true;
            print = true;
            printParDEM = true;

            outputTable.insert(basic_Info, printBasic);
            outputTable.insert(iB_Info, printIb);
            outputTable.insert(DEM_Info, printDEM);
            outputTable.insert(addModel_Info, printAddModel);
            outputTable.insert(parallelDEM_Info, printParallelDEM);
        }
        // Destructors -----------------------------------------------------
        ~InfoHFDIBDEM(){}

        void setOutput(bool basicB, bool ibB, bool demB, bool addModelB, bool parDEMB)
        {
            outputTable[basic_Info] = basicB;
            outputTable[iB_Info] = ibB;
            outputTable[DEM_Info] = demB;
            outputTable[addModel_Info] = addModelB;
            outputTable[parallelDEM_Info] = parDEMB;
        }

        InfoHFDIBDEM& operator<<(const outputType& data)
        {
            if(outputTable[data]  && data != parallelDEM_Info )
            {
                print = true;
                printParDEM = false;
            }
            else if (outputTable[data] && data == parallelDEM_Info)
            {
                printParDEM = true;
                print = false;
            }
            else
            {
                print = false;
                printParDEM = false;
            }
            return *this;
        }

        template<typename T>
        InfoHFDIBDEM& operator<<(const T& data)
        {
            if(printParDEM && !print)
            {
                Pout << data;
            }
            if(print && !printParDEM)
            {
                Info << data;
            }
            return *this;
        }

        typedef ostream& (manip)(ostream&);
        InfoHFDIBDEM& operator<<(manip m)
        {
            if(print && !printParDEM)
            {
                Info << m;
            }

            if(printParDEM && !print)
            {
                Pout << m;
            }
            return *this;
        }
    };
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "declareExternVars.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
