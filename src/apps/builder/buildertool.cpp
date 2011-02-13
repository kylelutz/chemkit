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

#include "buildertool.h"

// --- Construction and Destruction ---------------------------------------- //
BuilderTool::BuilderTool(BuilderWindow *builder)
    : chemkit::GraphicsTool()
{
    m_builder = builder;
}

BuilderTool::~BuilderTool()
{
}

// --- Properties ---------------------------------------------------------- //
BuilderWindow* BuilderTool::builder() const
{
    return m_builder;
}

chemkit::MoleculeEditor* BuilderTool::editor() const
{
    return builder()->editor();
}

QWidget* BuilderTool::settingsWidget()
{
    return 0;
}

// --- Event Handlers ------------------------------------------------------ //
void BuilderTool::cut()
{
}

void BuilderTool::copy()
{
}

void BuilderTool::paste()
{
}

void BuilderTool::del()
{
}

// --- Protected Methods --------------------------------------------------- //
void BuilderTool::setCanCut(bool canCut)
{
    builder()->setCanCut(canCut);
}

void BuilderTool::setCanCopy(bool canCopy)
{
    builder()->setCanCopy(canCopy);
}

void BuilderTool::setCanPaste(bool canPaste)
{
    builder()->setCanPaste(canPaste);
}

void BuilderTool::setCanDelete(bool canDelete)
{
    builder()->setCanDelete(canDelete);
}
