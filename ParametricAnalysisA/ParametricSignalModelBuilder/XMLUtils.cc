#include "XMLUtils.h"
#include <TList.h>
#include <TClass.h>
#include <TXMLAttr.h>

#include <iostream>
#include <stdlib.h>

#include <stdexcept>

#include <boost/regex.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>


using std::cout;
using std::cerr;
using std::endl;
//----------------------------------------------------------------------

std::vector<TXMLNode *>
XMLUtils::getChildren(TXMLNode *parent, const std::string &childName)
{
  std::vector<TXMLNode *> matchingChildren;
  for (TXMLNode *child = parent->GetChildren(); child != NULL; child = child->GetNextNode())
  {
    if (child->GetNodeType() != TXMLNode::kXMLElementNode)
      continue;

    if (childName != child->GetNodeName())
      // name does not match
      continue;

    matchingChildren.push_back(child);


  } // loop over all children

  return matchingChildren;
}

//----------------------------------------------------------------------

TXMLNode*
XMLUtils::getSingleChild(TXMLNode *parent, const std::string &childName)
{
  std::vector<TXMLNode *> matchingChildren = getChildren(parent, childName);
  if (matchingChildren.size() != 1)
  {
    cerr << "expected exactly one child named '" << childName << "' of parent '" << parent->GetNodeName() << "' but found " << matchingChildren.size() << endl;
    exit(1);
  }

  // cout << "found child " << matchingChildren[0]->GetNodeName() << endl;
  return matchingChildren[0];
}

//----------------------------------------------------------------------

TXMLNode*
XMLUtils::getOptionalSingleChild(TXMLNode *parent, const std::string &childName)
{
  std::vector<TXMLNode *> matchingChildren = getChildren(parent, childName);
  if (matchingChildren.size() > 1)
  {
    cerr << "expected at most one child named '" << childName << "' of parent '" << parent->GetNodeName() << "' but found " << matchingChildren.size() << endl;
    exit(1);
  }

  if (matchingChildren.size() == 0)
    return NULL;

  return matchingChildren[0];
}

//----------------------------------------------------------------------

TXMLAttr*
XMLUtils::getAttributeNode(TXMLNode *node, const std::string &attrName)
{
  TList *attributes = node->GetAttributes();
  if (attributes == NULL)
    return NULL;

  for (int i = 0; i < attributes->GetSize(); ++i)
    {
      TXMLAttr *attrNode = (TXMLAttr*) attributes->At(i);

      if (attrNode->GetName() == attrName)
        return attrNode;
    } // loop over all attributes

  return NULL;
}

//----------------------------------------------------------------------

std::string
XMLUtils::getAttribute(TXMLNode *node, const std::string &attrName)
{
  TXMLAttr *attrNode = getAttributeNode(node, attrName);

  if (attrNode != NULL)
    return attrNode->GetValue();

  cerr << "attribute '" << attrName << "' not found, exiting" << endl;
  exit(1);
}

//----------------------------------------------------------------------

std::string
XMLUtils::getAttribute(TXMLNode *node, const std::string &attrName, const std::string &defaultValue)
{
  TXMLAttr *attrNode = getAttributeNode(node, attrName);

  if (attrNode != NULL)
    return attrNode->GetValue();
  else
    return defaultValue;
}

//----------------------------------------------------------------------
std::vector<TXMLNode *>
XMLUtils::getElementChildren(TXMLNode *parent)
{
  std::vector<TXMLNode *> matchingChildren;
  for (TXMLNode *child = parent->GetChildren(); child != NULL; child = child->GetNextNode())
  {
    if (child->GetNodeType() != TXMLNode::kXMLElementNode)
      continue;

    matchingChildren.push_back(child);

  } // loop over all children

  return matchingChildren;
}

//----------------------------------------------------------------------
std::string
XMLUtils::getTextContent(TXMLNode *node)
{
  const char *raw_text = node->GetText();
  if (raw_text == NULL)
    throw std::invalid_argument(std::string("no text found in node ") + node->GetNodeName());

  return raw_text;
}

//----------------------------------------------------------------------

#include <boost/algorithm/string.hpp>

bool
XMLUtils::getBooleanContent(TXMLNode *node)
{
  std::string text = getTextContent(node);
  boost::trim(text);
  boost::to_lower(text);

  if (text == "true")
    return true;
  if (text == "false")
    return false;

  throw std::invalid_argument(std::string("expected either 'true' or 'false' in node ") + node->GetNodeName());
}

//----------------------------------------------------------------------

bool
XMLUtils::getBooleanContent(TXMLNode *parentNode, const std::string &childNodeName, bool defaultValue)
{
  std::vector<TXMLNode *> children = getChildren(parentNode, childNodeName);
  if (children.size() < 1)
    return defaultValue;

  if (children.size() > 1)
  {
    cerr << "more than one child with name '" << childNodeName << "' found " << endl;
    exit(1);
  }

  return getBooleanContent(children[0]);
}
//----------------------------------------------------------------------
#include <boost/lexical_cast.hpp>

int
XMLUtils::getIntegerContent(TXMLNode *node)
{
  // TODO: we should insist that there are no children other than
  // text nodes (and comment nodes)
  std::string text = getTextContent(node);
  boost::trim(text); // is this necessary before boost::lexical_cast ?

  return boost::lexical_cast<int>(text);
}
//----------------------------------------------------------------------

double
XMLUtils::getDoubleContent(TXMLNode *parentNode, const std::string &childNodeName, double defaultValue)
{
  std::vector<TXMLNode *> children = getChildren(parentNode, childNodeName);
  if (children.size() < 1)
    return defaultValue;

  if (children.size() > 1)
  {
    cerr << "more than one child with name '" << childNodeName << "' found " << endl;
    exit(1);
  }

  std::string text = getTextContent(parentNode);
  boost::trim(text);
  return boost::lexical_cast<double>(text);

}

//----------------------------------------------------------------------

std::string
XMLUtils::getStringContent(TXMLNode *parentNode, const std::string &childNodeName, const std::string &defaultValue)
{
  std::vector<TXMLNode *> children = getChildren(parentNode, childNodeName);
  if (children.size() < 1)
    return defaultValue;

  if (children.size() > 1)
  {
    cerr << "more than one child with name '" << childNodeName << "' found " << endl;
    exit(1);
  }

  std::string text = getTextContent(children[0]);
  boost::trim(text);
  return text;
}


//----------------------------------------------------------------------
std::set<unsigned> 
XMLUtils::getUnsignedSet(TXMLNode *parentNode, const std::string &childNodeName, const std::set<unsigned> &defaultValue)
{
  std::vector<TXMLNode *> children = getChildren(parentNode, childNodeName);
  if (children.size() < 1)
    return defaultValue;

  if (children.size() > 1)
  {
    cerr << "more than one child with name '" << childNodeName << "' found " << endl;
    exit(1);
  }

  // get the text
  std::string text = getTextContent(children[0]);
  boost::trim(text);

  // split the text (should be a comma separated list of integers)
  std::vector <std::string> parts;
  boost::algorithm::split_regex(parts, text,
                                boost::regex( "\\s*,\\s*"));

  std::set<unsigned> retval;
  BOOST_FOREACH(std::string part, parts)
  {
    std::cout << "PART='" << part << "'" << std::endl;
    retval.insert(boost::lexical_cast<unsigned>(part));
  }

  return retval;  
}
