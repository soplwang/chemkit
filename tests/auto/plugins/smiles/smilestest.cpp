/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
**
** This file is part of chemkit. For more information see
** <http://www.chemkit.org>.
**
** chemkit is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** chemkit is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with chemkit. If not, see <http://www.gnu.org/licenses/>.
**
******************************************************************************/

#include <QtTest>

#include <chemkit/molecule.h>
#include <chemkit/lineformat.h>
#include <chemkit/chemicalfile.h>
#include <chemkit/chemicalfileformat.h>

const QString dataPath = "../../../data/";

class SmilesTest : public QObject
{
    Q_OBJECT

    private slots:
        void initTestCase();

        // molecule tests
        void acenaphthylene();
        void aceticAcid();
        void adenine();
        void alanine();
        void ampicillin();
        void anthracene();
        void anthraquinone();
        void arsabenzene();
        void arsole();
        void aspirin();
        void aziridine();
        void azulene();
        void benzene();
        void benzofuran();
        void benzofurazan();
        void benzyne();
        void binol();
        void biphenyl();
        void biphenylene();
        void biperiden();
        void borinine();
        void borole();
        void buckminsterfullerene();
        void butene();
        void caffeine();
        void camphor();
        void carbazole();
        void cholesterol();
        void chrysene();
        void cinnoline();
        void colchicine();
        void copperSulfate();
        void corannulene();
        void coronene();
        void cubane();
        void cyanide();
        void cytosine();
        void decalin();
        void dibenzofuran();
        void dichloroethene();
        void dihydrogen();
        void dinitrogen();
        void ethane();
        void fluorenone();
        void folate();
        void furan();
        void furazan();
        void glucose();
        void guanine();
        void heavyWater();
        void histidine();
        void hydride();
        void hydronium();
        void ibuprofen();
        void indazole();
        void indene();
        void indole();
        void indolizine();
        void ipratropium();
        void isobutane();
        void isoindene();
        void isoindole();
        void melatonin();
        void naphthalene();
        void nicotine();
        void nitrobenzene();
        void ovalene();
        void oxazole();
        void pentacene();
        void pentalene();
        void perylene();
        void phenanthrene();
        void phenothiazine();
        void phenoxazine();
        void phosphole();
        void phosphorine();
        void phthalimide();
        void porphin();
        void proline();
        void proton();
        void purine();
        void pyranium();
        void pyrazole();
        void pyrene();
        void pyridazine();
        void pyridine();
        void pyrimidine();
        void pyrrole();
        void quinoline();
        void rhodizonicAcid();
        void selenophene();
        void sodiumChloride();
        void stilbene();
        void sulfurHexafluoride();
        void taxol();
        void tetralin();
        void tetraphenylene();
        void thiamin();
        void thiirane();
        void thiophene();
        void thujone();
        void thymine();
        void triazole();
        void triphenylene();
        void tropone();
        void tryptophan();
        void uracil();
        void vanillin();

        // feature tests
        void addHydrogens();
        void isotope();
        void kekulize();
        void quadrupleBond();

        // invalid tests
        void extraParenthesis();
        void invalidAtom();
        void wildcardAtom();

        // file tests
        void herg();
        void cox2();

    private:
        void COMPARE_SMILES(const chemkit::Molecule *molecule, const QString &smiles);
};

void SmilesTest::initTestCase()
{
    QVERIFY(chemkit::LineFormat::formats().contains("smiles"));
    QVERIFY(chemkit::ChemicalFileFormat::formats().contains("smi"));
}

void SmilesTest::COMPARE_SMILES(const chemkit::Molecule *molecule, const QString &smiles)
{
    chemkit::Molecule moleculeFromSmiles(smiles, "smiles");

    bool equal = molecule->equals(&moleculeFromSmiles, chemkit::Molecule::CompareAromaticity);
    if(!equal){
        qDebug() << "Actual SMILES: " << molecule->formula("smiles");
        qDebug() << "Actual formula: " << moleculeFromSmiles.formula();
        qDebug() << "Expected SMILES: " << smiles;
        qDebug() << "Expected formula: " << molecule->formula();
    }
    QVERIFY(equal);
}

// --- Molecule Tests ------------------------------------------------------ //
void SmilesTest::acenaphthylene()
{
    chemkit::Molecule molecule("c3cc1cccc2C=Cc(c12)c3", "smiles");
    QCOMPARE(molecule.formula(), QString("C12H8"));
    QCOMPARE(molecule.ringCount(), 3);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::aceticAcid()
{
    chemkit::Molecule molecule("CC(=O)O", "smiles");
    QCOMPARE(molecule.formula(), QString("C2H4O2"));
    QCOMPARE(molecule.bondCount(), 7);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::adenine()
{
    chemkit::Molecule molecule("n1c(c2c(nc1)ncn2)N", "smiles");
    QCOMPARE(molecule.formula(), QString("C5H5N5"));
    QCOMPARE(molecule.ringCount(), 2);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::alanine()
{
    chemkit::Molecule molecule("O=C(O)[C@H](N)C", "smiles");
    QCOMPARE(molecule.formula(), QString("C3H7NO2"));
    QCOMPARE(molecule.bondCount(), 12);

    foreach(const chemkit::Atom *atom, molecule.atoms()){
        if(atom->is(chemkit::Atom::Carbon) && atom->isBondedTo(chemkit::Atom::Nitrogen)){
            QVERIFY(atom->chirality() == chemkit::Atom::R);
        }
    }

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::ampicillin()
{
    chemkit::Molecule molecule("O=C(O)[C@@H]2N3C(=O)[C@@H](NC(=O)[C@@H]"
                               "(c1ccccc1)N)[C@H]3SC2(C)C", "smiles");
    QCOMPARE(molecule.formula(), QString("C16H19N3O4S"));
    QCOMPARE(molecule.ringCount(), 3);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::anthracene()
{
    chemkit::Molecule molecule("c1ccc2cc3ccccc3cc2c1", "smiles");
    QCOMPARE(molecule.formula(), QString("C14H10"));
    QCOMPARE(molecule.ringCount(), 3);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::anthraquinone()
{
    chemkit::Molecule molecule("O=C2c1ccccc1C(=O)c3ccccc23", "smiles");
    QCOMPARE(molecule.formula(), QString("C14H8O2"));
    QCOMPARE(molecule.ringCount(), 3);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::arsabenzene()
{
    chemkit::Molecule molecule("[as]1ccccc1", "smiles");
    QCOMPARE(molecule.formula(), QString("C5H5As"));
    QCOMPARE(molecule.ringCount(), 1);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::arsole()
{
    chemkit::Molecule molecule("c1[as]ccc1", "smiles");
    QCOMPARE(molecule.formula(), QString("C4H5As"));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::aspirin()
{
    chemkit::Molecule molecule("O=C(Oc1ccccc1C(=O)O)C", "smiles");
    QCOMPARE(molecule.formula(), QString("C9H8O4"));
    QCOMPARE(molecule.ringCount(), 1);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::aziridine()
{
    chemkit::Molecule molecule("N1CC1", "smiles");
    QCOMPARE(molecule.formula(), QString("C2H5N"));
    QCOMPARE(molecule.ringCount(), 1);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
    COMPARE_SMILES(&molecule, "C1NC1");
}

void SmilesTest::azulene()
{
    chemkit::Molecule molecule("c1cccc2cccc2c1", "smiles");
    QCOMPARE(molecule.formula(), QString("C10H8"));
    QCOMPARE(molecule.ringCount(), 2);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::benzene()
{
    chemkit::Molecule molecule("c1ccccc1", "smiles");
    QCOMPARE(molecule.formula(), QString("C6H6"));
    QCOMPARE(molecule.ringCount(), 1);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
    COMPARE_SMILES(&molecule, "C1=CC=CC=C1");
}

void SmilesTest::benzofuran()
{
    chemkit::Molecule molecule("o2c1ccccc1cc2", "smiles");
    QCOMPARE(molecule.formula(), QString("C8H6O"));
    QCOMPARE(molecule.ringCount(), 2);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::benzofurazan()
{
    chemkit::Molecule molecule("n1onc2ccccc12", "smiles");
    QCOMPARE(molecule.formula(), QString("C6H4N2O"));
    QCOMPARE(molecule.ringCount(), 2);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::benzyne()
{
    chemkit::Molecule molecule("C\\1#C\\C=C/C=C/1", "smiles");
    QCOMPARE(molecule.formula(), QString("C6H4"));
    QCOMPARE(molecule.ringCount(), 1);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::binol()
{
    chemkit::Molecule molecule("Oc1c(c2c(O)ccc3c2cccc3)c(cccc4)c4cc1", "smiles");
    QCOMPARE(molecule.formula(), QString("C20H14O2"));
    QCOMPARE(molecule.ringCount(), 4);

    foreach(const chemkit::Ring *ring, molecule.rings()){
        QCOMPARE(ring->size(), 6);
        QCOMPARE(ring->isAromatic(), true);
    }

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::biphenyl()
{
    chemkit::Molecule molecule("c1ccccc1(c2ccccc2)", "smiles");
    QCOMPARE(molecule.formula(), QString("C12H10"));
    QCOMPARE(molecule.ringCount(), 2);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::biphenylene()
{
    chemkit::Molecule molecule("c3cc2c1c(cccc1)c2cc3", "smiles");
    QCOMPARE(molecule.formula(), QString("C12H8"));
    QCOMPARE(molecule.ringCount(), 3);

    foreach(const chemkit::Ring *ring, molecule.rings()){
        if(ring->size() == 6){
            QCOMPARE(ring->isAromatic(), true);
        }
        else if(ring->size() == 4){
            QCOMPARE(ring->isAromatic(), false);
        }
    }

    //COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::biperiden()
{
    chemkit::Molecule molecule("OC(c1ccccc1)(CCN2CCCCC2)C4C3\\C=C/C(C3)C4", "smiles");
    QCOMPARE(molecule.formula(), QString("C21H29NO"));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::borinine()
{
    chemkit::Molecule molecule("b1ccccc1", "smiles");
    QCOMPARE(molecule.formula(), QString("C5H5B"));
    QCOMPARE(molecule.ringCount(), 1);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::borole()
{
    chemkit::Molecule molecule("C1=CC=CB1", "smiles");
    QCOMPARE(molecule.formula(), QString("C4H5B"));
    QCOMPARE(molecule.ringCount(), 1);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::buckminsterfullerene()
{
    chemkit::Molecule molecule("c12c3c4c5c1c6c7c8c2c9c%10c3c%11c%12c4c%13c%14"
                               "c5c%15c6c%16c7c%17c%18c8c9c%19c%20c%10c%11c%21"
                               "c%22c%12c%13c%23c%24c%14c%15c%25c%16c%26c%17"
                               "c%27c%18c%19c%28c%20c%21c%29c%22c%23c%30c%24"
                               "c%25c%26c%31c%27c%28c%29c%30%31", "smiles");
    QCOMPARE(molecule.formula(), QString("C60"));
}

void SmilesTest::butene()
{
    // cis butene
    chemkit::Molecule cis("C(=C\\C)\\C", "smiles");
    QCOMPARE(cis.formula(), QString("C4H8"));

    // trans butene
    chemkit::Molecule trans("C(=C/C)\\C", "smiles");
    QCOMPARE(trans.formula(), QString("C4H8"));
}

void SmilesTest::caffeine()
{
    chemkit::Molecule molecule("O=C2N(c1ncn(c1C(=O)N2C)C)C", "smiles");
    QCOMPARE(molecule.formula(), QString("C8H10N4O2"));
    QCOMPARE(molecule.ringCount(), 2);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::camphor()
{
    chemkit::Molecule molecule("O=C1CC2CCC1(C)C2(C)C", "smiles");
    QCOMPARE(molecule.formula(), QString("C10H16O"));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::carbazole()
{
    chemkit::Molecule molecule("c1cccc3c1c2c(cccc2)n3", "smiles");
    QCOMPARE(molecule.ringCount(), 3);

    foreach(const chemkit::Ring *ring, molecule.rings()){
        if(ring->contains(chemkit::Atom::Nitrogen)){
            QCOMPARE(ring->size(), 5);
        }
        else{
            QCOMPARE(ring->size(), 6);
        }

        QCOMPARE(ring->isAromatic(), true);
    }

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::cholesterol()
{
    chemkit::Molecule molecule("O[C@@H]4C/C3=C/C[C@@H]1[C@H](CC[C@]2([C@H]1CC"
                               "[C@@H]2[C@H](C)CCCC(C)C)C)[C@@]3(C)CC4", "smiles");
    QCOMPARE(molecule.formula(), QString("C27H46O"));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::chrysene()
{
    chemkit::Molecule molecule("c4c1c(ccc2ccccc12)c3ccccc3c4", "smiles");
    QCOMPARE(molecule.formula(), QString("C18H12"));
    QCOMPARE(molecule.ringCount(), 4);

    foreach(const chemkit::Ring *ring, molecule.rings()){
        QCOMPARE(ring->size(), 6);
        QCOMPARE(ring->isAromatic(), true);
    }

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::cinnoline()
{
    chemkit::Molecule molecule("n1nccc2ccccc12", "smiles");
    QCOMPARE(molecule.formula(), QString("C8H6N2"));
    QCOMPARE(molecule.ringCount(), 2);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::colchicine()
{
    chemkit::Molecule molecule("O=C(N[C@@H]3C\\1=C\\C(=O)C(\\OC)=C/C=C/1c2c"
                               "(cc(OC)c(OC)c2OC)CC3)C", "smiles");
    QCOMPARE(molecule.formula(), QString("C22H25NO6"));
    QCOMPARE(molecule.ringCount(), 3);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::copperSulfate()
{
    chemkit::Molecule molecule("[Cu+2].[O-]S(=O)(=O)[O-]", "smiles");
    QCOMPARE(molecule.formula(), QString("CuO4S"));
    QCOMPARE(molecule.bondCount(), 4);
    QCOMPARE(molecule.fragmentCount(), 2);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::corannulene()
{
    chemkit::Molecule molecule("c16ccc2ccc3ccc5c4c(c1c2c34)c(cc5)cc6", "smiles");
    QCOMPARE(molecule.formula(), QString("C20H10"));
    QCOMPARE(molecule.ringCount(), 6);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::coronene()
{
    chemkit::Molecule molecule("c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67", "smiles");
    QCOMPARE(molecule.formula(), QString("C24H12"));
    QCOMPARE(molecule.ringCount(), 7);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::cubane()
{
    chemkit::Molecule molecule("C12C3C4C1C5C2C3C45", "smiles");
    QCOMPARE(molecule.formula(), QString("C8H8"));
    QCOMPARE(molecule.bondCount(), 20);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::cyanide()
{
    chemkit::Molecule molecule("C#N", "smiles");
    QCOMPARE(molecule.formula(), QString("CHN"));
    QCOMPARE(molecule.bondCount(), 2);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::cytosine()
{
    chemkit::Molecule molecule("O=C1/N=C\\C=C(\\N)N1", "smiles");
    QCOMPARE(molecule.formula(), QString("C4H5N3O"));
    QCOMPARE(molecule.ringCount(), 1);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::decalin()
{
    chemkit::Molecule molecule("C1CCC2CCCCC2C1", "smiles");
    QCOMPARE(molecule.formula(), QString("C10H18"));
    QCOMPARE(molecule.ringCount(), 2);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::dibenzofuran()
{
    chemkit::Molecule molecule("o2c1ccccc1c3c2cccc3", "smiles");
    QCOMPARE(molecule.formula(), QString("C12H8O"));
    QCOMPARE(molecule.ringCount(), 3);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);
    QCOMPARE(molecule.rings()[2]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::dichloroethene()
{
    chemkit::Molecule molecule("Cl[C@H]=CCl", "smiles");
    QCOMPARE(molecule.formula(), QString("C2H2Cl2"));
    QCOMPARE(molecule.bondCount(), 5);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::dihydrogen()
{
    chemkit::Molecule molecule("[H][H]", "smiles");
    QCOMPARE(molecule.formula(), QString("H2"));
    QCOMPARE(molecule.bondCount(), 1);

    QCOMPARE(molecule.formula("smiles"), QString("[H][H]"));
    //COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::dinitrogen()
{
    chemkit::Molecule molecule("N#N", "smiles");
    QCOMPARE(molecule.formula(), QString("N2"));

    QCOMPARE(molecule.formula("smiles"), QString("N#N"));
    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::ethane()
{
    chemkit::Molecule molecule("CC", "smiles");
    QCOMPARE(molecule.formula(), QString("C2H6"));

    QCOMPARE(molecule.formula("smiles"), QString("CC"));
    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::fluorenone()
{
    chemkit::Molecule molecule("O=C3c1ccccc1c2c3cccc2", "smiles");
    QCOMPARE(molecule.formula(), QString("C13H8O"));
    QCOMPARE(molecule.ringCount(), 3);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::folate()
{
    chemkit::Molecule molecule("O=C(O)[C@@H](NC(=O)c1ccc(cc1)NCc2nc3c"
                               "(nc2)N/C(=N\\C3=O)N)CCC(=O)O", "smiles");
    QCOMPARE(molecule.formula(), QString("C19H19N7O6"));
    QCOMPARE(molecule.ringCount(), 3);

    //COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::furan()
{
    chemkit::Molecule molecule("o1cccc1", "smiles");
    QCOMPARE(molecule.formula(), QString("C4H4O"));
    QCOMPARE(molecule.ringCount(), 1);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::furazan()
{
    chemkit::Molecule molecule("n1oncc1", "smiles");
    QCOMPARE(molecule.formula(), QString("C2H2N2O"));
    QCOMPARE(molecule.ringCount(), 1);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::glucose()
{
    chemkit::Molecule molecule("OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)1", "smiles");
    QCOMPARE(molecule.formula(), QString("C6H12O6"));
    QCOMPARE(molecule.ringCount(), 1);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::guanine()
{
    chemkit::Molecule molecule("NC1=Nc2[nH]cnc2C(=O)N1", "smiles");
    QCOMPARE(molecule.formula(), QString("C5H5N5O"));
    QCOMPARE(molecule.ringCount(), 2);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::heavyWater()
{
    chemkit::Molecule molecule("[2H]O[2H]", "smiles");
    QCOMPARE(molecule.formula(), QString("H2O"));
    QCOMPARE(molecule.bondCount(), 2);

    foreach(const chemkit::Atom *atom, molecule.atoms()){
        if(atom->is(chemkit::Atom::Hydrogen)){
            QCOMPARE(atom->massNumber(), 2);
        }
    }

    QCOMPARE(molecule.formula("smiles"), QString("[2H]O[2H]"));
    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::histidine()
{
    chemkit::Molecule molecule("N[C@@H](Cc1[nH]cnc1)C(O)=O", "smiles");
    QCOMPARE(molecule.formula(), QString("C6H9N3O2"));
    QCOMPARE(molecule.ringCount(), 1);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::hydride()
{
    chemkit::Molecule molecule("[H-]", "smiles");
    QCOMPARE(molecule.formula(), QString("H"));
    QCOMPARE(molecule.bondCount(), 0);
    //QCOMPARE(molecule.atom(0)->formalCharge(), -1);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::hydronium()
{
    chemkit::Molecule molecule("[OH3+]", "smiles");
    QCOMPARE(molecule.formula(), QString("H3O"));
    QCOMPARE(molecule.bondCount(), 3);

    foreach(const chemkit::Atom *atom, molecule.atoms()){
        if(atom->is(chemkit::Atom::Oxygen)){
            QCOMPARE(atom->formalCharge(), 1);
        }
    }

    QCOMPARE(molecule.formula("smiles"), QString("[OH3+]"));
    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::ibuprofen()
{
    chemkit::Molecule molecule("CC(C(=O)O)c1ccc(CC(C)C)cc1", "smiles");
    QCOMPARE(molecule.formula(), QString("C13H18O2"));
    QCOMPARE(molecule.ringCount(), 1);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::indazole()
{
    chemkit::Molecule molecule("n2cc1ccccc1n2", "smiles");
    QCOMPARE(molecule.formula(), QString("C7H6N2"));
    QCOMPARE(molecule.ringCount(), 2);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::indene()
{
    chemkit::Molecule molecule("c1cccc2c1\\C=C/C2", "smiles");
    QCOMPARE(molecule.formula(), QString("C9H8"));
    QCOMPARE(molecule.ringCount(), 2);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::indole()
{
    chemkit::Molecule molecule("c1cccc2c1ccn2", "smiles");
    QCOMPARE(molecule.formula(), QString("C8H7N"));
    QCOMPARE(molecule.ringCount(), 2);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::indolizine()
{
    chemkit::Molecule molecule("c1ccc2ccccn12", "smiles");
    QCOMPARE(molecule.formula(), QString("C8H7N"));
    QCOMPARE(molecule.ringCount(), 2);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::ipratropium()
{
    chemkit::Molecule molecule("O=C(OC2CC1[N+](C)(C(CC1)C2)C(C)C)C(c3ccccc3)CO", "smiles");
    QCOMPARE(molecule.formula(), QString("C20H30NO3"));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::isobutane()
{
    chemkit::Molecule molecule("CC(C)C", "smiles");
    QCOMPARE(molecule.formula(), QString("C4H10"));
    QCOMPARE(molecule.bondCount(), 13);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::isoindene()
{
    chemkit::Molecule molecule("C12=CCC=C1C=CC=C2", "smiles");
    QCOMPARE(molecule.formula(), QString("C9H8"));
    QCOMPARE(molecule.ringCount(), 2);

    foreach(const chemkit::Ring *ring, molecule.rings()){
        if(ring->size() == 5){
            QCOMPARE(ring->isAromatic(), false);
        }
        else if(ring->size() == 6){
            QCOMPARE(ring->isAromatic(), true);
        }
    }

    //COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::isoindole()
{
    chemkit::Molecule molecule("c1cccc2c1cnc2", "smiles");
    QCOMPARE(molecule.formula(), QString("C8H7N"));
    QCOMPARE(molecule.ringCount(), 2);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::melatonin()
{
    chemkit::Molecule molecule("O=C(NCCc2c1cc(OC)ccc1nc2)C", "smiles");
    QCOMPARE(molecule.formula(), QString("C13H16N2O2"));
    QCOMPARE(molecule.ringCount(), 2);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::naphthalene()
{
    chemkit::Molecule molecule("c1ccc2ccccc2c1", "smiles");
    QCOMPARE(molecule.formula(), QString("C10H8"));
    QCOMPARE(molecule.ringCount(), 2);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::nicotine()
{
    chemkit::Molecule molecule("CN1CCC[C@H]1c2cccnc2", "smiles");
    QCOMPARE(molecule.formula(), QString("C10H14N2"));
    QCOMPARE(molecule.ringCount(), 2);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::nitrobenzene()
{
    chemkit::Molecule molecule("[O-][N+](=O)c1ccccc1", "smiles");
    QCOMPARE(molecule.formula(), QString("C6H5NO2"));
    QCOMPARE(molecule.ringCount(), 1);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::ovalene()
{
    chemkit::Molecule molecule("c1cc2c3c4c1ccc5cc6c7c8c(ccc9=c8c1c(cc9)cc"
                               "(c3c1c7c54)cc2)cc6", "smiles");
    QCOMPARE(molecule.formula(), QString("C32H14"));
    QCOMPARE(molecule.ringCount(), 10);
}

void SmilesTest::oxazole()
{
    chemkit::Molecule molecule("n1ccoc1", "smiles");
    QCOMPARE(molecule.formula(), QString("C3H3NO"));
    QCOMPARE(molecule.ringCount(), 1);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::pentacene()
{
    chemkit::Molecule molecule("c45cc3cc2cc1ccccc1cc2cc3cc4cccc5", "smiles");
    QCOMPARE(molecule.formula(), QString("C22H14"));
    QCOMPARE(molecule.ringCount(), 5);

    foreach(const chemkit::Ring *ring, molecule.rings()){
        QCOMPARE(ring->size(), 6);
        QCOMPARE(ring->isAromatic(), true);
    }

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::pentalene()
{
    chemkit::Molecule molecule("c1cc2cccc2c1", "smiles");
    QCOMPARE(molecule.formula(), QString("C8H6"));
    QCOMPARE(molecule.ringCount(), 2);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::perylene()
{
    chemkit::Molecule molecule("c1ccc5cccc4c5c1c2cccc3cccc4c23", "smiles");
    QCOMPARE(molecule.formula(), QString("C20H12"));
    QCOMPARE(molecule.ringCount(), 5);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::phenanthrene()
{
    chemkit::Molecule molecule("c1ccc2c(c1)ccc3ccccc32", "smiles");
    QCOMPARE(molecule.formula(), QString("C14H10"));
    QCOMPARE(molecule.ringCount(), 3);

    foreach(const chemkit::Ring *ring, molecule.rings()){
        QCOMPARE(ring->size(), 6);
        QCOMPARE(ring->isAromatic(), true);
    }

    //COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::phenothiazine()
{
    chemkit::Molecule molecule("c1ccc2Sc3ccccc3Nc2c1", "smiles");
    QCOMPARE(molecule.formula(), QString("C12H9NS"));
    QCOMPARE(molecule.ringCount(), 3);

    //COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::phenoxazine()
{
    chemkit::Molecule molecule("O2c1ccccc1Nc3c2cccc3", "smiles");
    QCOMPARE(molecule.formula(), QString("C12H9NO"));
    QCOMPARE(molecule.ringCount(), 3);

    //COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::phosphole()
{
    chemkit::Molecule molecule("c1cccp1", "smiles");
    QCOMPARE(molecule.formula(), QString("C4H5P"));
    QCOMPARE(molecule.ringCount(), 1);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::phosphorine()
{
    chemkit::Molecule molecule("p1ccccc1", "smiles");
    QCOMPARE(molecule.formula(), QString("C5H5P"));
    QCOMPARE(molecule.ringCount(), 1);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::phthalimide()
{
    chemkit::Molecule molecule("O=C2c1ccccc1C(=O)N2", "smiles");
    QCOMPARE(molecule.formula(), QString("C8H5NO2"));
    QCOMPARE(molecule.ringCount(), 2);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::porphin()
{
    chemkit::Molecule molecule("c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3", "smiles");
    QCOMPARE(molecule.formula(), QString("C20H14N4"));
    QCOMPARE(molecule.ringCount(), 5);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::proline()
{
    chemkit::Molecule molecule("O=C(O)C1NCCC1", "smiles");
    QCOMPARE(molecule.formula(), QString("C5H9NO2"));
    QCOMPARE(molecule.ringCount(), 1);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::proton()
{
    chemkit::Molecule molecule("[H+]", "smiles");
    QCOMPARE(molecule.formula(), QString("H"));
    QCOMPARE(molecule.bondCount(), 0);

    QCOMPARE(molecule.formula("smiles"), QString("[H+]"));
    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::purine()
{
    chemkit::Molecule molecule("n1cc2c(nc1)ncn2", "smiles");
    QCOMPARE(molecule.formula(), QString("C5H4N4"));
    QCOMPARE(molecule.ringCount(), 2);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::pyranium()
{
    chemkit::Molecule molecule("[o+]1ccccc1", "smiles");
    QCOMPARE(molecule.formula(), QString("C5H5O"));
    QCOMPARE(molecule.ringCount(), 1);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    foreach(const chemkit::Atom *atom, molecule.atoms()){
        if(atom->is(chemkit::Atom::Oxygen)){
            QCOMPARE(atom->formalCharge(), 1);
        }
    }

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::pyrazole()
{
    chemkit::Molecule molecule("n1cccn1", "smiles");
    QCOMPARE(molecule.formula(), QString("C3H4N2"));
    QCOMPARE(molecule.ringCount(), 1);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::pyrene()
{
    chemkit::Molecule molecule("c3ccc2ccc1cccc4c1c2c3cc4", "smiles");
    QCOMPARE(molecule.formula(), QString("C16H10"));
    QCOMPARE(molecule.ringCount(), 4);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::pyridazine()
{
    chemkit::Molecule molecule("n1ncccc1", "smiles");
    QCOMPARE(molecule.formula(), QString("C4H4N2"));
    QCOMPARE(molecule.ringCount(), 1);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::pyridine()
{
    chemkit::Molecule molecule("n1ccccc1", "smiles");
    QCOMPARE(molecule.formula(), QString("C5H5N"));
    QCOMPARE(molecule.ringCount(), 1);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::pyrimidine()
{
    chemkit::Molecule molecule("n1cccnc1", "smiles");
    QCOMPARE(molecule.formula(), QString("C4H4N2"));
    QCOMPARE(molecule.ringCount(), 1);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::pyrrole()
{
    chemkit::Molecule molecule("n1cccc1", "smiles");
    QCOMPARE(molecule.formula(), QString("C4H5N"));
    QCOMPARE(molecule.ringCount(), 1);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::quinoline()
{
    chemkit::Molecule molecule("n1cccc2ccccc12", "smiles");
    QCOMPARE(molecule.formula(), QString("C9H7N"));
    QCOMPARE(molecule.ringCount(), 2);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::rhodizonicAcid()
{
    chemkit::Molecule molecule("O=C1C(/O)=C(/O)C(=O)C(=O)C1=O", "smiles");
    QCOMPARE(molecule.formula(), QString("C6H2O6"));
    QCOMPARE(molecule.ringCount(), 1);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::selenophene()
{
    chemkit::Molecule molecule("[se]1cccc1", "smiles");
    QCOMPARE(molecule.formula(), QString("C4H4Se"));
    QCOMPARE(molecule.ringCount(), 1);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::sodiumChloride()
{
    chemkit::Molecule molecule("[Na+].[Cl-]", "smiles");
    QCOMPARE(molecule.formula(), QString("ClNa"));
    QCOMPARE(molecule.bondCount(), 0);
    QCOMPARE(molecule.fragmentCount(), 2);

    QCOMPARE(molecule.formula("smiles"), QString("[Na+].[Cl-]"));
    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::stilbene()
{
    chemkit::Molecule molecule("c2(\\C=C\\c1ccccc1)ccccc2", "smiles");
    QCOMPARE(molecule.formula(), QString("C14H12"));
    QCOMPARE(molecule.bondCount(), 27);
    QCOMPARE(molecule.ringCount(), 2);

    foreach(const chemkit::Ring *ring, molecule.rings()){
        QCOMPARE(ring->size(), 6);
        QCOMPARE(ring->isAromatic(), true);
    }

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::sulfurHexafluoride()
{
    chemkit::Molecule molecule("FS(F)(F)(F)(F)F", "smiles");
    QCOMPARE(molecule.formula(), QString("F6S"));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::taxol()
{
    chemkit::Molecule molecule("O=C(N[C@@H](c1ccccc1)[C@@H](O)C(=O)O[C@H]5C"
                               "[C@@]6(O)[C@@H](OC(=O)c2ccccc2)[C@H]3[C@@](C)"
                               "([C@@H](O)C[C@H]4OC[C@@]34OC(C)=O)C(=O)[C@H]"
                               "(OC(C)=O)\\C(=C5/C)[C@]6(C)C)c7ccccc7", "smiles");
    QCOMPARE(molecule.formula(), QString("C47H51NO14"));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::tetraphenylene()
{
    chemkit::Molecule molecule("c5cc4c1c(cccc1)c2ccccc2c3ccccc3c4cc5", "smiles");
    QCOMPARE(molecule.formula(), QString("C24H16"));
    QCOMPARE(molecule.ringCount(), 5);
}

void SmilesTest::tetralin()
{
    chemkit::Molecule molecule("c1ccc2c(c1)CCCC2", "smiles");
    QCOMPARE(molecule.formula(), QString("C10H12"));
    QCOMPARE(molecule.ringCount(), 2);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::thiamin()
{
    chemkit::Molecule molecule("n1c(c(cnc1C)C[n+]2c(c(sc2)CCO)C)N", "smiles");
    QCOMPARE(molecule.formula(), QString("C12H16N4OS"));
    QCOMPARE(molecule.ringCount(), 2);

    //COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::thiirane()
{
    chemkit::Molecule molecule("C1CS1", "smiles");
    QCOMPARE(molecule.formula(), QString("C2H4S"));
    QCOMPARE(molecule.ringCount(), 1);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::thiophene()
{
    chemkit::Molecule molecule("s1cccc1", "smiles");
    QCOMPARE(molecule.formula(), QString("C4H4S"));
    QCOMPARE(molecule.ringCount(), 1);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::thujone()
{
    chemkit::Molecule molecule("C[C@@H]([C@@H](C2)[C@]2([C@@H](C)C)C1)C1=O", "smiles");
    QCOMPARE(molecule.formula(), QString("C10H16O"));
    QCOMPARE(molecule.ringCount(), 2);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::thymine()
{
    chemkit::Molecule molecule("O=C1\\C(=C/NC(=O)N1)C", "smiles");
    QCOMPARE(molecule.formula(), QString("C5H6N2O2"));
    QCOMPARE(molecule.ringCount(), 1);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::triazole()
{
    chemkit::Molecule molecule("n1ccnn1", "smiles");
    QCOMPARE(molecule.formula(), QString("C2H3N3"));
    QCOMPARE(molecule.ringCount(), 1);
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::triphenylene()
{
    chemkit::Molecule molecule("c4cc3c1c(cccc1)c2ccccc2c3cc4", "smiles");
    QCOMPARE(molecule.formula(), QString("C18H12"));
    QCOMPARE(molecule.ringCount(), 4);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
    COMPARE_SMILES(&molecule, "C1(C=CC=C3)=C3C(C=CC=C4)=C4C2=C1C=CC=C2");
}

void SmilesTest::tropone()
{
    chemkit::Molecule molecule("C1=CC=CC(=O)C=C1", "smiles");
    QCOMPARE(molecule.formula(), QString("C7H6O"));
    QCOMPARE(molecule.ringCount(), 1);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::tryptophan()
{
    chemkit::Molecule molecule("N[C@@H](Cc1c2ccccc2nc1)C(O)=O", "smiles");
    QCOMPARE(molecule.formula(), QString("C11H12N2O2"));
    QCOMPARE(molecule.ringCount(), 2);

    foreach(const chemkit::Ring *ring, molecule.rings()){
        if(ring->contains(chemkit::Atom::Nitrogen)){
            QCOMPARE(ring->size(), 5);
        }
        else{
            QCOMPARE(ring->size(), 6);
        }

        QCOMPARE(ring->isAromatic(), true);
    }

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::uracil()
{
    chemkit::Molecule molecule("O=C1\\C=C/NC(=O)N1", "smiles");
    QCOMPARE(molecule.formula(), QString("C4H4N2O2"));
    QCOMPARE(molecule.ringCount(), 1);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::vanillin()
{
    chemkit::Molecule molecule("O=CC1=CC(OC)=C(O)C=C1", "smiles");
    QCOMPARE(molecule.formula(), QString("C8H8O3"));
    QCOMPARE(molecule.ringCount(), 1);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

// --- Feature Tests ------------------------------------------------------- //
void SmilesTest::addHydrogens()
{
    chemkit::LineFormat *format = chemkit::LineFormat::create("smiles");
    QVERIFY(format);

    // default is true
    QCOMPARE(format->option("add-hydrogens").toBool(), true);
    chemkit::Molecule *molecule = format->read("C");
    QCOMPARE(molecule->formula(), QString("CH4"));
    delete molecule;

    format->setOption("add-hydrogens", false);
    molecule = format->read("C");
    QCOMPARE(molecule->formula(), QString("C"));
    delete molecule;

    delete format;
}

void SmilesTest::isotope()
{
    chemkit::LineFormat *format = chemkit::LineFormat::create("smiles");
    QVERIFY(format);

    chemkit::Molecule *molecule = format->read("[14CH4]");
    QVERIFY(molecule);
    QCOMPARE(molecule->formula(), QString("CH4"));

    foreach(const chemkit::Atom *atom, molecule->atoms()){
        if(atom->is(chemkit::Atom::Carbon)){
            QCOMPARE(atom->massNumber(), 14);
        }
    }

    delete molecule;

    molecule = format->read("[238U]");
    QVERIFY(molecule);
    QCOMPARE(molecule->formula(), QString("U"));
    QCOMPARE(molecule->atom(0)->massNumber(), 238);

    delete molecule;
    delete format;
}

void SmilesTest::kekulize()
{
    chemkit::LineFormat *format = chemkit::LineFormat::create("smiles");
    QVERIFY(format);

    // default is false
    QCOMPARE(format->option("kekulize").toBool(), false);

    chemkit::Molecule benzene("c1ccccc1", "smiles");
    QCOMPARE(format->write(&benzene), QString("c1ccccc1"));

    format->setOption("kekulize", true);
    QCOMPARE(format->write(&benzene), QString("C1=CC=CC=C1"));

    delete format;
}

void SmilesTest::quadrupleBond()
{
    chemkit::LineFormat *format = chemkit::LineFormat::create("smiles");
    QVERIFY(format);

    chemkit::Molecule *molecule = format->read("C$C");
    QVERIFY(molecule);
    QCOMPARE(molecule->formula(), QString("C2"));
    QCOMPARE(molecule->bondCount(), 1);
    QCOMPARE(molecule->bonds()[0]->order(), 4);

    delete molecule;
    delete format;
}

// --- Invalid Tests ------------------------------------------------------- //
void SmilesTest::extraParenthesis()
{
    chemkit::LineFormat *format = chemkit::LineFormat::create("smiles");
    QVERIFY(format);

    chemkit::Molecule *molecule = format->read("C(C=O))C");
    QVERIFY(molecule == 0);
    QVERIFY(format->errorString().isEmpty() == false);

    delete format;
}

void SmilesTest::invalidAtom()
{
    chemkit::LineFormat *format = chemkit::LineFormat::create("smiles");
    QVERIFY(format);

    chemkit::Molecule *molecule = format->read("CCX");
    QVERIFY(molecule == 0);
    QVERIFY(format->errorString().isEmpty() == false);

    delete format;
}

void SmilesTest::wildcardAtom()
{
    chemkit::LineFormat *format = chemkit::LineFormat::create("smiles");
    QVERIFY(format);

    chemkit::Molecule *molecule = format->read("C*C");
    QVERIFY(molecule == 0);
    QVERIFY(format->errorString().isEmpty() == false);

    delete format;
}

// --- File Tests ---------------------------------------------------------- //
void SmilesTest::herg()
{
    chemkit::ChemicalFile file(dataPath + "herg.smi");

    bool ok = file.read();
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString();
    QVERIFY(ok);

    QCOMPARE(file.moleculeCount(), 31);
    QCOMPARE(file.molecule(0)->name(), QString("Amitriptyline"));
    QCOMPARE(file.molecule(0)->formula(), QString("C20H23N"));
    QCOMPARE(file.molecule(30)->name(), QString("Verapamil"));
    QCOMPARE(file.molecule(30)->formula(), QString("C27H38N2O4"));
}

void SmilesTest::cox2()
{
    chemkit::ChemicalFile file(dataPath + "cox2.smi");

    bool ok = file.read();
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString();
    QVERIFY(ok);

    QCOMPARE(file.moleculeCount(), 128);
    QCOMPARE(file.molecule(0)->formula(), QString("C13H18N2O5S"));
    QCOMPARE(file.molecule(2)->formula(), QString("C16H13F2NO3S2"));
    QCOMPARE(file.molecule(127)->formula(), QString("C21H19NO5S"));
}

QTEST_APPLESS_MAIN(SmilesTest)
#include "smilestest.moc"