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

#include "buildtool.h"

#include <chemkit/vector.h>
#include <chemkit/element.h>
#include <chemkit/periodictabledialog.h>

// --- Construction and Destruction ---------------------------------------- //
BuildTool::BuildTool(BuilderWindow *builder)
    : QObject(),
      BuilderTool(builder)
{
    m_element = chemkit::Atom::Carbon;
    m_bondOrder = chemkit::Bond::Single;
    m_adjustHydrogens = true;

    m_intialAtom = 0;
    m_movingAtom = 0;
    m_bondingAtom = 0;
    m_newBond = 0;

    // elements to display in element selector
    m_elements.append(chemkit::Atom::Carbon);
    m_elements.append(chemkit::Atom::Nitrogen);
    m_elements.append(chemkit::Atom::Oxygen);
    m_elements.append(chemkit::Atom::Chlorine);
    m_elements.append(chemkit::Atom::Bromine);
    m_elements.append(chemkit::Atom::Hydrogen);
    m_elements.append(chemkit::Atom::Phosphorus);
    m_elements.append(chemkit::Atom::Sulfur);
}

BuildTool::~BuildTool()
{
}

/// Sets the current element to \p element.
void BuildTool::setElement(const chemkit::Element &element)
{
    if(element.isValid()){
        m_element = element;
    }

    if(m_elementSelector){
        if(m_elements.contains(element.atomicNumber())){
            m_elementSelector->setCurrentIndex(m_elements.indexOf(element.atomicNumber()));
        }
        else if(m_addedElements.contains(element.atomicNumber())){
            m_elementSelector->setCurrentIndex(m_elementSelector->findText(element.name()));
        }
        else{
            m_elementSelector->removeItem(m_elementSelector->count() - 1);
            m_elementSelector->addItem(element.name(), element.atomicNumber());
            m_elementSelector->addItem("Other", -1);
            m_elementSelector->update();
            m_addedElements.append(element.atomicNumber());

            m_elementSelector->setCurrentIndex(m_elementSelector->count() - 2);
        }
    }
}

/// Returns the current element.
chemkit::Element BuildTool::element() const
{
    return m_element;
}

/// Sets the current bond order to \p bondOrder.
void BuildTool::setBondOrder(int bondOrder)
{
    m_bondOrder = bondOrder;
}

/// Returns the current bond order.
int BuildTool::bondOrder() const
{
    return m_bondOrder;
}

// --- Settings ------------------------------------------------------------ //
QWidget* BuildTool::settingsWidget()
{
    QWidget *widget = new QWidget;

    QFormLayout *layout = new QFormLayout;

    // element selector
    m_elementSelector = new QComboBox;
    foreach(int element, m_elements){
        m_elementSelector->addItem(chemkit::Element(element).name(), element);
    }
    if(!m_addedElements.isEmpty()){
        foreach(int element, m_addedElements){
            m_elementSelector->addItem(chemkit::Element(element).name(), element);
        }
    }
    m_elementSelector->addItem("Other", -1);
    connect(m_elementSelector, SIGNAL(currentIndexChanged(int)), SLOT(elementSelectorChanged(int)));
    layout->addRow("Element:", m_elementSelector);

    // bond order selector
    m_bondOrderSelector = new QComboBox;
    m_bondOrderSelector->addItem("Single");
    m_bondOrderSelector->addItem("Double");
    m_bondOrderSelector->addItem("Triple");
    m_bondOrderSelector->setCurrentIndex(m_bondOrder - 1);
    connect(m_bondOrderSelector, SIGNAL(currentIndexChanged(int)), SLOT(bondOrderSelectorChanged(int)));
    layout->addRow("Bond Order:", m_bondOrderSelector);

    // add hydrogens checkbox
    m_addHydrogensCheckBox = new QCheckBox("Auto Add Hydrogens");
    m_addHydrogensCheckBox->setChecked(m_adjustHydrogens);
    connect(m_addHydrogensCheckBox, SIGNAL(stateChanged(int)), SLOT(addHydrogensChanged(int)));
    layout->addRow(m_addHydrogensCheckBox);

    widget->setLayout(layout);

    return widget;
}

// --- Events -------------------------------------------------------------- //
void BuildTool::mousePressEvent(QMouseEvent *event)
{
    beginMoleculeEdit();

    // left button
    if(event->button() == Qt::LeftButton){
        chemkit::GraphicsItem *item = view()->itemAt(event->x(), event->y());

        // item under cursor
        if(item){
            // atom under cursor
            if(item->type() == chemkit::GraphicsItem::AtomItem){
                chemkit::GraphicsAtomItem *atomItem = static_cast<chemkit::GraphicsAtomItem *>(item);
                chemkit::Atom *atom = const_cast<chemkit::Atom *>(atomItem->atom());
                m_intialElement = atom->atomicNumber();

                if(atom->atomicNumber() != m_element.atomicNumber()){
                    setAtomAtomicNumber(atom, m_element.atomicNumber());
                }

                m_intialAtom = atom;
            }

            // bond under cursor
            else if(item->type() == chemkit::GraphicsItem::BondItem){
                chemkit::GraphicsBondItem *bondItem = static_cast<chemkit::GraphicsBondItem *>(item);
                chemkit::Bond *bond = const_cast<chemkit::Bond *>(bondItem->bond());
                setBondOrder(bond, (bond->order() % 3) + 1);
            }
        }

        // no item under cursor
        else{
            // add new atom
            chemkit::Atom *atom = addAtom(m_element.atomicNumber());
            chemkit::GraphicsPoint position = view()->unproject(event->x(), event->y(), editor()->molecule()->center());
            setAtomPosition(atom, position.toPoint());
            m_intialAtom = atom;
            m_intialElement = m_element.atomicNumber();
        }

        m_movingAtom = 0;
        m_bondingAtom = 0;
        m_newBond = 0;
    }

    // right button
    else if(event->button() == Qt::RightButton){
        chemkit::GraphicsItem *item = view()->itemAt(event->x(), event->y());

        if(item){
            if(item->type() == chemkit::GraphicsItem::AtomItem){
                chemkit::GraphicsAtomItem *atomItem = static_cast<chemkit::GraphicsAtomItem *>(item);
                chemkit::Atom *atom = const_cast<chemkit::Atom *>(atomItem->atom());
                removeAtom(atom);
            }
            else if(item->type() == chemkit::GraphicsItem::BondItem){
                chemkit::GraphicsBondItem *bondItem = static_cast<chemkit::GraphicsBondItem *>(item);
                chemkit::Bond *bond = const_cast<chemkit::Bond *>(bondItem->bond());
                removeBond(bond);
            }
        }
    }

    view()->update();
}

void BuildTool::mouseMoveEvent(QMouseEvent *event)
{
    if(!m_intialAtom){
        return;
    }

    chemkit::GraphicsItem *item = 0;

    if(m_movingAtom){
        QList<chemkit::GraphicsItem *> items = view()->itemsAt(event->x(), event->y());

        if(items.isEmpty()){
            item = 0;
        }
        else if(items.size() == 1){
            item = items[0];
        }
        else{
            chemkit::GraphicsItem *nearestItem = items[0];

            // if the nearest item that we intersect is the atom item for
            // the currently moving atom, then select the next item beneath it.
            if(nearestItem->type() == chemkit::GraphicsItem::AtomItem &&
               static_cast<chemkit::GraphicsAtomItem *>(nearestItem)->atom() == m_movingAtom){
                item = items[1];
            }
            else{
                item = items[0];
            }
        }
    }
    else{
        item = view()->itemAt(event->x(), event->y());
    }

    // cursor over nothing
    if(!item){
        if(!m_movingAtom){
            setAtomAtomicNumber(m_intialAtom, m_intialElement);
            m_movingAtom = addAtom(m_element.atomicNumber());
            addBond(m_intialAtom, m_movingAtom, bondOrder());
            chemkit::GraphicsPoint position = view()->unproject(event->x(), event->y(), m_intialAtom->position());
            setAtomPosition(m_movingAtom, position.toPoint());

            if(m_newBond){
                removeBond(m_newBond);
                m_newBond = 0;
                m_bondingAtom = 0;
            }
        }
        else{
            chemkit::GraphicsPoint newPosition = view()->unproject(event->x(), event->y(), m_movingAtom->position());
            setAtomPosition(m_movingAtom, newPosition.toPoint());
        }
    }

    // cursor over atom item
    else if(item->type() == chemkit::GraphicsItem::AtomItem){
        chemkit::GraphicsAtomItem *atomItem = static_cast<chemkit::GraphicsAtomItem *>(item);
        chemkit::Atom *atom = const_cast<chemkit::Atom *>(atomItem->atom());

        // over intial atom
        if(atom == m_intialAtom){
            if(m_movingAtom){
                removeAtom(m_movingAtom);
                m_movingAtom = 0;
                setAtomAtomicNumber(m_intialAtom, m_element.atomicNumber());
            }
        }
        // over moving atom
        else if(atom == m_movingAtom){
            chemkit::GraphicsPoint newPosition = view()->unproject(event->x(), event->y(), m_movingAtom->position());
            setAtomPosition(m_movingAtom, newPosition.toPoint());
        }
        // over new atom
        else{
            if(m_movingAtom){
                removeAtom(m_movingAtom);
                m_movingAtom = 0;
            }

            if(m_newBond && atom != m_bondingAtom){
                removeBond(m_newBond);
                m_bondingAtom = 0;
            }

            if(!m_intialAtom->isBondedTo(atom)){
                m_newBond = addBond(atom, m_intialAtom);
                setBondOrder(m_newBond, m_bondOrder);
                m_bondingAtom = atom;
            }
        }
    }

    view()->update();
}

void BuildTool::mouseReleaseEvent(QMouseEvent *event)
{
    if(event->button() == Qt::LeftButton){
        m_intialAtom = 0;
        m_movingAtom = 0;
        m_bondingAtom = 0;
        m_newBond = 0;
    }

    if(m_adjustHydrogens){
        foreach(chemkit::Atom *atom, m_modifiedAtoms){
            adjustHydrogens(atom);
        }
    }
    m_modifiedAtoms.clear();

    endMoleculeEdit();
}

// --- Slots --------------------------------------------------------------- //
void BuildTool::elementSelectorChanged(int index)
{
    int atomicNumber = m_elementSelector->itemData(index).toInt();
    if(atomicNumber == -1){
        const chemkit::Element &element = chemkit::PeriodicTableDialog::getElement(builder(), "Select Element");

        if(element.isValid()){
            setElement(element);
        }
        else{
            setElement(chemkit::Atom::Carbon);
        }
    }
    else{
        setElement(atomicNumber);
    }
}

void BuildTool::bondOrderSelectorChanged(int index)
{
    setBondOrder(index + 1);
}

void BuildTool::addHydrogensChanged(int state)
{
    if(state == Qt::Checked){
        m_adjustHydrogens = true;
    }
    else{
        m_adjustHydrogens = false;
    }
}

// --- Internal Methods ---------------------------------------------------- //
void BuildTool::beginMoleculeEdit()
{
    builder()->beginMoleculeEdit();
}

void BuildTool::endMoleculeEdit()
{
    // do hydrogen adjustment
    if(m_adjustHydrogens){
        foreach(chemkit::Atom *atom, m_modifiedAtoms){
            adjustHydrogens(atom);
        }
    }
    m_modifiedAtoms.clear();

    builder()->endMoleculeEdit();
}

chemkit::Atom* BuildTool::addAtom(int atomicNumber)
{
    chemkit::Atom *atom = editor()->addAtom(atomicNumber);
    m_modifiedAtoms.insert(atom);
    return atom;
}

void BuildTool::removeAtom(chemkit::Atom *atom)
{
    foreach(chemkit::Atom *neighbor, atom->neighbors()){
        m_modifiedAtoms.insert(neighbor);
    }

    m_modifiedAtoms.remove(atom);

    editor()->removeAtom(atom);
}

void BuildTool::setAtomAtomicNumber(chemkit::Atom *atom, int atomicNumber)
{
    m_modifiedAtoms.insert(atom);
    editor()->setAtomAtomicNumber(atom, atomicNumber);
}

void BuildTool::setAtomPosition(chemkit::Atom *atom, const chemkit::Point &position)
{
    editor()->setAtomPosition(atom, position);
}

chemkit::Bond* BuildTool::addBond(chemkit::Atom *a, chemkit::Atom *b, int order)
{
    m_modifiedAtoms.insert(a);
    m_modifiedAtoms.insert(b);
    return editor()->addBond(a, b, order);
}

void BuildTool::removeBond(chemkit::Bond *bond)
{
    m_modifiedAtoms.insert(bond->atom1());
    m_modifiedAtoms.insert(bond->atom2());
    editor()->removeBond(bond);
}

void BuildTool::setBondOrder(chemkit::Bond *bond, int order)
{
    m_modifiedAtoms.insert(bond->atom1());
    m_modifiedAtoms.insert(bond->atom2());
    editor()->setBondOrder(bond, order);
}

void BuildTool::adjustHydrogens(chemkit::Atom *atom)
{
    // remove lone hydrogens
    if(atom->is(chemkit::Atom::Hydrogen) && atom->neighborCount() < 2){
        editor()->removeAtom(atom);
    }

    // add hydrogens
    else if(atom->valence() < atom->expectedValence()){
        while(atom->valence() < atom->expectedValence()){
            chemkit::Atom *hydrogen = editor()->addAtom(chemkit::Atom::Hydrogen);
            editor()->setAtomPosition(hydrogen, atom->position().movedBy(chemkit::Vector::randomUnitVector()));
            editor()->addBond(atom, hydrogen);
        }
    }

    // remove hydrogens
    else if(atom->valence() > atom->expectedValence()){
        foreach(chemkit::Atom *neighbor, atom->neighbors()){
            if(neighbor->isTerminalHydrogen()){
                editor()->removeAtom(neighbor);
                m_modifiedAtoms.remove(neighbor);

                if(atom->valence() == atom->expectedValence())
                    break;
            }
        }
    }
}