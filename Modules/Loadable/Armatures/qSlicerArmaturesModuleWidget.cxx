/*=========================================================================

  Program: Bender

  Copyright (c) Kitware Inc.

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0.txt

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

=========================================================================*/

// Qt includes

// SlicerQt includes
#include "qSlicerArmaturesModuleWidget.h"
#include "ui_qSlicerArmaturesModule.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_Armatures
class qSlicerArmaturesModuleWidgetPrivate: public Ui_qSlicerArmaturesModule
{
public:
  qSlicerArmaturesModuleWidgetPrivate();
};

//-----------------------------------------------------------------------------
// qSlicerArmaturesModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerArmaturesModuleWidgetPrivate::qSlicerArmaturesModuleWidgetPrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerArmaturesModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerArmaturesModuleWidget::qSlicerArmaturesModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerArmaturesModuleWidgetPrivate )
{
}

//-----------------------------------------------------------------------------
qSlicerArmaturesModuleWidget::~qSlicerArmaturesModuleWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidget::setup()
{
  Q_D(qSlicerArmaturesModuleWidget);
  d->setupUi(this);
  this->Superclass::setup();
}

