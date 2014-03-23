#include <QObject>

class StepHook: public QObject
{
  Q_OBJECT
public:
  StepHook();
  ~StepHook(){}
public slots:
  void step();
public:
  class Internal;
  Internal* Pimpl;
private:
  Q_DISABLE_COPY(StepHook);
};
