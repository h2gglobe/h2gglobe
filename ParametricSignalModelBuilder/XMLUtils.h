#ifndef _XMLUtils_h
#define _XMLUtils_h

#include <TXMLNode.h>
#include <TXMLAttr.h>

#include <vector>
#include <string>

class XMLUtils
{
public:
  
  /** @return all children with the given name */
  static std::vector<TXMLNode *> getChildren(TXMLNode *parent, const std::string &childName);

  /*  exits the program / throws an exception if no or more than one matching
   *  child was found.
   */
  static TXMLNode* getSingleChild(TXMLNode *parent, const std::string &childName);

  /** @return the child if present, NULL if no such child was found and throws
   *  an exception / stops the program if more than one child was found.
   */
  static TXMLNode* getOptionalSingleChild(TXMLNode *parent, const std::string &childName);

  /** @return the attribute node or NULL if no such attribute exists */
  static TXMLAttr* getAttributeNode(TXMLNode *node, const std::string &attrName);
  
  /** exits/throws an exception if the attribute is not found */
  static std::string getAttribute(TXMLNode *node, const std::string &attrName);

  /** @return the value of an attribute or the given default value if the attribute does not exist */
  static std::string getAttribute(TXMLNode *node, const std::string &attrName, const std::string &defaultValue);
  
  /** @return all (direct) children of the given node which are element nodes
   *  (i.e. no attributes, comment nodes etc.)
   */
  static std::vector<TXMLNode *> getElementChildren(TXMLNode *parent);

  /** Calls node->GetText() and insists that there is a text or throws an exception otherwise
   */
  static std::string getTextContent(TXMLNode *node);

  //--------------------

  /** convenience method for reading configuration files: returns the boolean value of a note of the form:
   *
   *    <node>true</node>   or <node>false</node>
   *
   *  Throws an exception if the text in node is not 'true' nor 'false' (where the comparison
   *  is NOT case-sensitive)
   */
  static bool getBooleanContent(TXMLNode *node);

  /** similar to the other getBooleanContent(..) function but returns a default value if the given
   *  node is not present.
   */
  static bool getBooleanContent(TXMLNode *parentNode, const std::string &childNodeName, bool defaultValue);

  /** gets the text content of this node and requires that it corresponds to an integer and
   *  returns the corresponding value.
   */
  static int getIntegerContent(TXMLNode *node);

  /** gets the text content of this node and requires that it corresponds to an integer and
   *  returns the corresponding value.
   */
  static double getDoubleContent(TXMLNode *parentNode, const std::string &childNodeName, double defaultValue);
  //--------------------

};


#endif
