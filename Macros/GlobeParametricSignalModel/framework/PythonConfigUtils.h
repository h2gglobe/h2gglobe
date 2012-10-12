#ifndef _PythonConfigUtils_h
#define _PythonConfigUtils_h

#include "pythonFwd.h"


#include <string>
#include <vector>


class PythonConfigUtils
{
public:

  /** checks whether the given dict has a key with the given name. If yes, ensures that
   *  the value is also a dict, if not returns NULL.
   */
  static PyObject *getOptionalDict(PyObject *parentDict, const std::string &keyName);

  static PyObject* getDictValue(PyObject *dict, const std::string &keyName);

  /** tries to find the given key in the given dict and ensures that this is a string. */
  static std::string getDictStringValue(PyObject *dict, const std::string &keyName);

  /** similar to the other getDictStringValue method but returns the given default value
   *  if the key is not present in the dict.
   */
  static std::string getDictStringValue(PyObject *dict, const std::string &keyName, const std::string &defaultValue);

  /** tries to find the given key in the given dict and ensures that this is a boolean. */
  static bool getDictBooleanValue(PyObject *dict, const std::string &keyName);

  /** similar to the other getDictBooleanValue(..) method but returns the default value
   *  if the key is not found.
   */
  static bool getDictBooleanValue(PyObject *dict, const std::string &keyName, bool defaultValue);

  /** tries to find the given key in the given dict and ensures that this is an integer. */
  static long getDictIntValue(PyObject *dict, const std::string &keyName);

  static long getDictIntValue(PyObject *dict, const std::string &keyName, long defaultValue);

  static double getDictDoubleValue(PyObject *dict, const std::string &keyName);

  /** tries to find the given key in the given dict and ensures that this is a double. */
  static double getDictDoubleValue(PyObject *dict, const std::string &keyName, double defaultValue);

  static std::vector<int> getDictIntListValue(PyObject *dict, const std::string &keyName, const std::vector<int> &defaultValue);

  static std::vector<std::string> getDictStringListValue(PyObject *dict, const std::string &keyName, const std::vector<std::string> &defaultValue);

  /** @return a pointer to a (newly allocated) double or NULL if this key
   *  is not found in the dict.
   */
  static double *getDictDoubleValuePointer(PyObject *dict, const std::string &keyName);

  static std::string *getDictStringValuePointer(PyObject *dict, const std::string &keyName);

  static bool *getDictBooleanValuePointer(PyObject *dict, const std::string &keyName);

  static double castNumberToFloat(PyObject *obj);

  static std::string castToString(PyObject *obj);
};


#endif
