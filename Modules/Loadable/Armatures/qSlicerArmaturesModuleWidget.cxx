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
#include <QDebug>
#include <QVector3d>

// Armatures includes
#include "qSlicerArmaturesModuleWidget.h"
#include "qSlicerArmaturesModuleWidget_p.h"
#include "ui_qSlicerArmaturesModule.h"
#include "vtkMRMLArmatureNode.h"
#include "vtkMRMLBoneNode.h"
#include "vtkSlicerArmaturesLogic.h"

// VTK includes
#include <vtkCollection.h>
#include <vtkNew.h>
#include <vtkStdString.h>

// Annotations includes
#include <qMRMLSceneAnnotationModel.h>
#include <vtkSlicerAnnotationModuleLogic.h>

// MRML includes
#include <vtkMRMLInteractionNode.h>
#include <vtkMRMLSelectionNode.h>

//-----------------------------------------------------------------------------
// qSlicerArmaturesModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerArmaturesModuleWidgetPrivate
::qSlicerArmaturesModuleWidgetPrivate(qSlicerArmaturesModuleWidget& object)
  : q_ptr(&object)
{
  this->ArmatureNode = 0;
  this->BoneNode = 0;
}

//-----------------------------------------------------------------------------
vtkSlicerArmaturesLogic* qSlicerArmaturesModuleWidgetPrivate::logic() const
{
  Q_Q(const qSlicerArmaturesModuleWidget);
  return vtkSlicerArmaturesLogic::SafeDownCast(q->logic());
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidgetPrivate
::setupUi(qSlicerWidget* armatureModuleWidget)
{
  Q_Q(qSlicerArmaturesModuleWidget);
  this->Superclass::setupUi(armatureModuleWidget);

  // Armatures
  QObject::connect(this->ArmatureNodeComboBox,
                   SIGNAL(currentNodeChanged(vtkMRMLNode*)),
                   q, SLOT(setMRMLArmatureNode(vtkMRMLNode*)));
  QObject::connect(this->ArmatureVisibilityCheckBox, SIGNAL(toggled(bool)),
                   q, SLOT(setArmatureVisibility(bool)));

  // Bones
  this->BonesTreeView->annotationModel()->setAnnotationsAreParent(true);
  this->BonesTreeView->setLogic(this->logic()->GetAnnotationsLogic());
  this->BonesTreeView->annotationModel()->setNameColumn(0);
  this->BonesTreeView->annotationModel()->setVisibilityColumn(0);

  this->BonesTreeView->annotationModel()->setCheckableColumn(-1);
  this->BonesTreeView->annotationModel()->setLockColumn(-1);
  this->BonesTreeView->annotationModel()->setEditColumn(-1);
  this->BonesTreeView->annotationModel()->setValueColumn(-1);
  this->BonesTreeView->annotationModel()->setTextColumn(-1);

  this->BonesTreeView->setHeaderHidden(true);

  QObject::connect(this->BonesTreeView,
                   SIGNAL(currentNodeChanged(vtkMRMLNode*)),
                   q, SLOT(setMRMLBoneNode(vtkMRMLNode*)));

  QAction* addBoneAction = new QAction("Add bone", this->BonesTreeView);
  this->BonesTreeView->prependNodeMenuAction(addBoneAction);
  this->BonesTreeView->prependSceneMenuAction(addBoneAction);
  QObject::connect(addBoneAction, SIGNAL(triggered()),
                   q, SLOT(addAndPlaceBone()));

  // Logic
  q->qvtkConnect(this->logic(), vtkCommand::ModifiedEvent,
                 q, SLOT(updateWidgetFromLogic()));

  QObject::connect(this->BonesTreeView,
    SIGNAL(currentNodeChanged(vtkMRMLNode*)),
    q, SLOT(onTreeNodeSelected(vtkMRMLNode*)));

  // -- Positions --
  QObject::connect(this->HeadCoordinatesWidget,
    SIGNAL(coordinatesChanged(double*)),
    q, SLOT(updateCurrentMRMLBoneNode()));
  QObject::connect(this->TailCoordinatesWidget,
    SIGNAL(coordinatesChanged(double*)),
    q, SLOT(updateCurrentMRMLBoneNode()));

  QObject::connect(this->BonePositionTypeComboBox,
    SIGNAL(currentIndexChanged(QString)),
    this, SLOT(onPositionTypeChanged()));
  QObject::connect(this->LengthSpinBox,
    SIGNAL(valueChanged(double)),
    this, SLOT(onDistanceChanged(double)));
  QObject::connect(this->DirectionCoordinatesWidget,
    SIGNAL(coordinatesChanged(double*)),
    this, SLOT(onDirectionChanged(double*)));

  // -- Display --
  this->BoneRepresentationComboBox->setCurrentIndex(2);

  QObject::connect(this->BoneNameLineEdit,
    SIGNAL(textChanged(const QString&)), q, SLOT(updateCurrentMRMLBoneNode()));
  QObject::connect(this->BoneRepresentationComboBox,
    SIGNAL(currentIndexChanged(int)), q, SLOT(updateCurrentMRMLBoneNode()));
  QObject::connect(this->BoneColorPickerButton,
    SIGNAL(colorChanged(QColor)), q, SLOT(updateCurrentMRMLBoneNode()));
  QObject::connect(this->BoneOpacitySlider,
    SIGNAL(valueChanged(double)), q, SLOT(updateCurrentMRMLBoneNode()));
  QObject::connect(this->BoneShowAxesComboBox,
    SIGNAL(currentIndexChanged(int)), q, SLOT(updateCurrentMRMLBoneNode()));
  QObject::connect(this->BoneShowParenthoodCheckBox,
    SIGNAL(stateChanged(int)), q, SLOT(updateCurrentMRMLBoneNode()));
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidgetPrivate::onPositionTypeChanged()
{
  bool oldHeadState = this->HeadCoordinatesWidget->blockSignals(true);
  bool oldTailState = this->TailCoordinatesWidget->blockSignals(true);

  this->setCoordinatesFromBoneNode(this->BoneNode);

  this->HeadCoordinatesWidget->blockSignals(oldHeadState);
  this->TailCoordinatesWidget->blockSignals(oldTailState);
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidgetPrivate::onDistanceChanged(double newDistance)
{
  QVector3D head(this->HeadCoordinatesWidget->coordinates()[0],
    this->HeadCoordinatesWidget->coordinates()[1],
    this->HeadCoordinatesWidget->coordinates()[2]);
  QVector3D direction(this->DirectionCoordinatesWidget->coordinates()[0],
    this->DirectionCoordinatesWidget->coordinates()[1],
    this->DirectionCoordinatesWidget->coordinates()[2]);

  QVector3D newTail = head + direction * newDistance;

  this->TailCoordinatesWidget->setCoordinates(
    newTail.x(), newTail.y(), newTail.z());
}

//-----------------------------------------------------------------------------
QVector3D qSlicerArmaturesModuleWidgetPrivate::direction()
{
  QVector3D head(this->HeadCoordinatesWidget->coordinates()[0],
    this->HeadCoordinatesWidget->coordinates()[1],
    this->HeadCoordinatesWidget->coordinates()[2]);
  QVector3D tail(this->TailCoordinatesWidget->coordinates()[0],
    this->TailCoordinatesWidget->coordinates()[1],
    this->TailCoordinatesWidget->coordinates()[2]);

  return QVector3D(tail - head).normalized();
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidgetPrivate::blockPositionsSignals(bool block)
{
  this->HeadCoordinatesWidget->blockSignals(block);
  this->TailCoordinatesWidget->blockSignals(block);
  this->LengthSpinBox->blockSignals(block);
  this->DirectionCoordinatesWidget->blockSignals(block);
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidgetPrivate::blockDisplaySignals(bool block)
{
  this->BoneNameLineEdit->blockSignals(block);
  this->BoneRepresentationComboBox->blockSignals(block);
  this->BoneColorPickerButton->blockSignals(block);
  this->BoneOpacitySlider->blockSignals(block);
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidgetPrivate
::onDirectionChanged(double* newDirection)
{
  QVector3D head(this->HeadCoordinatesWidget->coordinates()[0],
    this->HeadCoordinatesWidget->coordinates()[1],
    this->HeadCoordinatesWidget->coordinates()[2]);
  QVector3D direction(newDirection[0], newDirection[1], newDirection[2]);

  QVector3D newTail = head + direction * this->LengthSpinBox->value();

  this->TailCoordinatesWidget->setCoordinates(
    newTail.x(), newTail.y(), newTail.z());
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidgetPrivate
::updateArmatureWidget(vtkMRMLBoneNode* boneNode)
{
  this->updateHierarchy(boneNode);
  this->updatePositions(boneNode);
  this->updateDisplay(boneNode);
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidgetPrivate
::updateHierarchy(vtkMRMLBoneNode* boneNode)
{
  this->HierarchyCollapsibleGroupBox->setEnabled(boneNode != 0);
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidgetPrivate
::updatePositions(vtkMRMLBoneNode* boneNode)
{
  this->blockPositionsSignals(true);

  bool enableHeadWidget = false;
  bool enableTailWidget = false;

  this->setCoordinatesFromBoneNode(boneNode);

  QVector3D direction(this->direction());
  this->DirectionCoordinatesWidget->setCoordinates(
    direction.x(), direction.y(), direction.z());

  if (boneNode)
    {
    this->LengthSpinBox->setValue(boneNode->GetLength());

    if (boneNode->GetWidgetState() == vtkMRMLBoneNode::PlaceTail)
      {
      enableHeadWidget = true;
      }
    else if (boneNode->GetWidgetState() == vtkMRMLBoneNode::Rest)
      {
      enableHeadWidget = true;
      enableTailWidget = true;
      }
    }
  else
    {
    this->LengthSpinBox->setValue(0.0);
    }

  this->HeadCoordinatesWidget->setEnabled(enableHeadWidget);
  this->TailCoordinatesWidget->setEnabled(enableTailWidget);
  this->LengthSpinBox->setEnabled(enableTailWidget);
  this->DirectionCoordinatesWidget->setEnabled(enableTailWidget);

  this->updateAdvancedPositions(boneNode);

  this->PositionsCollapsibleGroupBox->setEnabled(boneNode != 0);
  this->blockPositionsSignals(false);
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidgetPrivate
::updateAdvancedPositions(vtkMRMLBoneNode* boneNode)
{
  Q_UNUSED(boneNode);
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidgetPrivate
::updateDisplay(vtkMRMLBoneNode* boneNode)
{
  if (boneNode)
    {
    this->blockDisplaySignals(true);
    this->BoneNameLineEdit->setText(QString(boneNode->GetBoneName()));

    this->BoneRepresentationComboBox->setCurrentIndex(
      boneNode->GetBoneRepresentationType());

    int rgb[3];
    boneNode->GetBoneColor(rgb);
    this->BoneColorPickerButton->setColor(
      QColor::fromRgb(rgb[0], rgb[1], rgb[2]));

    this->BoneOpacitySlider->setValue(boneNode->GetOpacity());
    this->blockDisplaySignals(false);
    }

  this->updateAdvancedDisplay(boneNode);

  this->BoneDisplayCollapsibleGroupBox->setEnabled(boneNode != 0);
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidgetPrivate
::updateAdvancedDisplay(vtkMRMLBoneNode* boneNode)
{
  if (boneNode)
    {
    this->BoneShowAxesComboBox->setCurrentIndex(
      boneNode->GetShowAxes());
    this->BoneShowParenthoodCheckBox->setChecked(
      boneNode->GetShowParenthood());
    }
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidgetPrivate
::setCoordinatesToBoneNode(vtkMRMLBoneNode* boneNode)
{
  if (!boneNode)
    {
    return;
    }

  if (boneNode->GetWidgetState() == vtkMRMLBoneNode::Rest)
    {
    if (this->BonePositionTypeComboBox->currentText() == "Local")
      {
      boneNode->SetLocalHeadRest(this->HeadCoordinatesWidget->coordinates());
      boneNode->SetLocalTailRest(this->TailCoordinatesWidget->coordinates());
      }
    else
      {
      boneNode->SetWorldHeadRest(this->HeadCoordinatesWidget->coordinates());
      boneNode->SetWorldTailRest(this->TailCoordinatesWidget->coordinates());
      }
    }
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidgetPrivate
::setCoordinatesFromBoneNode(vtkMRMLBoneNode* boneNode)
{
  double head[3] = {0.0, 0.0, 0.0};
  double tail[3] = {0.0, 0.0, 0.0};

  if (boneNode && boneNode->GetWidgetState() == vtkMRMLBoneNode::Rest)
    {
    if (this->BonePositionTypeComboBox->currentText() == "Local")
      {
      boneNode->GetLocalHeadRest(head);
      boneNode->GetLocalTailRest(tail);
      }
    else
      {
      boneNode->GetWorldHeadRest(head);
      boneNode->GetWorldTailRest(tail);
      }
    }
  else if (boneNode && boneNode->GetWidgetState() == vtkMRMLBoneNode::Pose)
    {
    if (this->BonePositionTypeComboBox->currentText() == "Local")
      {
      boneNode->GetLocalHeadPose(head);
      boneNode->GetLocalTailPose(tail);
      }
    else
      {
      boneNode->GetWorldHeadPose(head);
      boneNode->GetWorldTailPose(tail);
      }
    }

  this->HeadCoordinatesWidget->setCoordinates(head);
  this->TailCoordinatesWidget->setCoordinates(tail);
}

//-----------------------------------------------------------------------------
// qSlicerArmaturesModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerArmaturesModuleWidget::qSlicerArmaturesModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerArmaturesModuleWidgetPrivate(*this) )
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

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidget
::setMRMLArmatureNode(vtkMRMLArmatureNode* armatureNode)
{
  Q_D(qSlicerArmaturesModuleWidget);
  this->qvtkReconnect(d->ArmatureNode, armatureNode,
    vtkCommand::ModifiedEvent, this, SLOT(updateWidgetFromArmatureNode()));
  d->ArmatureNode = armatureNode;

  d->logic()->SetActiveArmature(armatureNode);
  this->onTreeNodeSelected(armatureNode);
  this->updateWidgetFromArmatureNode();
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidget
::setMRMLArmatureNode(vtkMRMLNode* armatureNode)
{
  this->setMRMLArmatureNode(vtkMRMLArmatureNode::SafeDownCast(armatureNode));
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidget::setArmatureVisibility(bool visible)
{
  Q_D(qSlicerArmaturesModuleWidget);
  vtkMRMLArmatureNode* armatureNode = this->mrmlArmatureNode();

  if (!armatureNode)
    {
    return;
    }

  armatureNode->SetVisibility(visible);

  vtkNew<vtkCollection> bones;
  armatureNode->GetAllBones(bones.GetPointer());
  for (int i = 0; i < bones->GetNumberOfItems(); ++i)
    {
    vtkMRMLBoneNode* boneNode
      = vtkMRMLBoneNode::SafeDownCast(bones->GetItemAsObject(i));

    if (boneNode)
      {
      boneNode->SetVisible(visible);
      }
    }
}

//-----------------------------------------------------------------------------
vtkMRMLArmatureNode* qSlicerArmaturesModuleWidget::mrmlArmatureNode()const
{
  Q_D(const qSlicerArmaturesModuleWidget);
  return vtkMRMLArmatureNode::SafeDownCast(
    d->ArmatureNodeComboBox->currentNode());
}

//-----------------------------------------------------------------------------
vtkMRMLBoneNode* qSlicerArmaturesModuleWidget::mrmlBoneNode()const
{
  Q_D(const qSlicerArmaturesModuleWidget);
  return vtkMRMLBoneNode::SafeDownCast(
    d->BonesTreeView->currentNode());
}

/*
//-----------------------------------------------------------------------------
vtkMRMLArmatureDisplayNode* qSlicerArmaturesModuleWidget
::mrmlArmatureDisplayNode()
{
  vtkMRMLArmatureNode* armatureNode = this->mrmlArmatureNode();
  return armatureNode ? armatureNode->GetArmatureDisplayNode() : 0;
}
*/

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidget
::setMRMLBoneNode(vtkMRMLBoneNode* boneNode)
{
  Q_D(qSlicerArmaturesModuleWidget);
  //d->logic()->SetActiveBone(boneNode);
  //if (boneNode == 0)
  //  {
  //  d->logic()->SetActiveArmature(this->mrmlArmatureNode());
  //  }
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidget
::setMRMLBoneNode(vtkMRMLNode* boneNode)
{
  Q_D(qSlicerArmaturesModuleWidget);
  this->setMRMLBoneNode(vtkMRMLBoneNode::SafeDownCast(boneNode));
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidget::addAndPlaceBone()
{
  Q_D(qSlicerArmaturesModuleWidget);
  vtkMRMLSelectionNode* selectionNode = vtkMRMLSelectionNode::SafeDownCast(
    this->mrmlScene()->GetNodeByID("vtkMRMLSelectionNodeSingleton"));
  vtkMRMLInteractionNode* interactionNode =
    vtkMRMLInteractionNode::SafeDownCast(
      this->mrmlScene()->GetNodeByID("vtkMRMLInteractionNodeSingleton"));
  if (!selectionNode || !interactionNode)
    {
    qCritical() << "Invalid scene, no interaction or selection node";
    return;
    }
  selectionNode->SetReferenceActiveAnnotationID("vtkMRMLBoneNode");
  interactionNode->SwitchToSinglePlaceMode();
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidget::updateWidgetFromLogic()
{
  Q_D(qSlicerArmaturesModuleWidget);
  vtkMRMLNode* activeNode = d->logic()->GetActiveBone();
  if (activeNode == 0)
    {
    activeNode = d->logic()->GetActiveArmature();
    }
  d->BonesTreeView->setCurrentNode(activeNode);
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidget::updateWidgetFromArmatureNode()
{
  Q_D(qSlicerArmaturesModuleWidget);

  if (!d->ArmatureNode)
    {
    return;
    }

  d->ArmatureVisibilityCheckBox->setChecked(d->ArmatureNode->GetVisibility());
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidget::updateWidgetFromBoneNode()
{
  Q_D(qSlicerArmaturesModuleWidget);
  d->updateArmatureWidget(d->BoneNode);
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidget::onTreeNodeSelected(vtkMRMLNode* node)
{
  Q_D(qSlicerArmaturesModuleWidget);

  vtkMRMLBoneNode* boneNode = vtkMRMLBoneNode::SafeDownCast(node);
  if (boneNode)
    {
    this->qvtkReconnect(d->BoneNode, boneNode, vtkCommand::ModifiedEvent,
      this, SLOT(updateWidgetFromBoneNode()));
    d->BoneNode = boneNode;
    }

  d->updateArmatureWidget(boneNode);
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesModuleWidget::updateCurrentMRMLBoneNode()
{
  Q_D(qSlicerArmaturesModuleWidget);

  if (!d->BoneNode)
    {
    return;
    }

  int wasModifying = d->BoneNode->StartModify();

  d->setCoordinatesToBoneNode(d->BoneNode);

  d->BoneNode->SetBoneName(d->BoneNameLineEdit->text().toStdString().c_str());

  d->BoneNode->SetBoneRepresentationType(
    d->BoneRepresentationComboBox->currentIndex());
  int rgb[3];
  rgb[0] = d->BoneColorPickerButton->color().red();
  rgb[1] = d->BoneColorPickerButton->color().green();
  rgb[2] = d->BoneColorPickerButton->color().blue();
  //qDebug()<<"SetBoneColor: "<<rgb[0]<<" "<<rgb[1]<<" "<<rgb[2];
  d->BoneNode->SetBoneColor(rgb);

  d->BoneNode->SetOpacity(d->BoneOpacitySlider->value());

  d->BoneNode->SetShowAxes(d->BoneShowAxesComboBox->currentIndex());

  d->BoneNode->SetShowParenthood(
    d->BoneShowParenthoodCheckBox->isChecked());

  d->BoneNode->EndModify(wasModifying);
}
