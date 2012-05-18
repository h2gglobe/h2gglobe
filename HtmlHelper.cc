#include "HtmlHelper.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>

#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"

//---------------------------------------------------------------------------------------------------------------------
template <typename T>
T& dereference(T* ptr) { return *ptr; }

//---------------------------------------------------------------------------------------------------------------------
template <typename T>
void delete_ptr(T * ptr) { if( ptr != 0 ) { delete ptr; } } 

//---------------------------------------------------------------------------------------------------------------------
#define RENDER(CLASS,TAG)                  	\
void CLASS::render(std::ostream & out )  {	\
	out << indent() << "<"  TAG << attributes_ << ">" \
	    << text_.str() << "\n";   		\
	out << objects_;			\
	out << unindent() << "</" TAG ">";	\
}

//---------------------------------------------------------------------------------------------------------------------
HtmlObject& HtmlObject::set(const std::string & n, const std::string & a)
{
	attributes_.push_back(std::make_pair(n,a));
	return *this;
}

//---------------------------------------------------------------------------------------------------------------------
class printAtt 
{
public:
	printAtt( std::ostream & out ) : out_(out) {}
	
	void operator() (const std::pair<std::string,std::string>& obj ) {
		out_ << " " << obj.first << "=\"" << obj.second << "\"";
	}

private:
	std::ostream & out_;
};
	
//---------------------------------------------------------------------------------------------------------------------
std::ostream & operator << ( std::ostream & out, const std::list<std::pair<std::string,std::string> > & obj)
{
	for_each(obj.begin(),obj.end(),printAtt(out));
	return out << ( obj.empty() ? "" : " " );
}

//---------------------------------------------------------------------------------------------------------------------
int HtmlObject::indent_(0);

//---------------------------------------------------------------------------------------------------------------------
std::string HtmlObject::indent()
{
	std::string ret;
	for(int ii=0; ii<indent_; ii++ ) { ret+= " "; }
	++indent_;
	return ret;
}

//---------------------------------------------------------------------------------------------------------------------
std::string HtmlObject::unindent()
{
	std::string ret;
	--indent_;
	for(int ii=0; ii<indent_; ii++ ) { ret+= " "; }
	return ret;
}

//---------------------------------------------------------------------------------------------------------------------
HtmlObject::HtmlObject()
	: own_(true)
{}

//---------------------------------------------------------------------------------------------------------------------
HtmlObject::~HtmlObject()
{
	if( own_ ) {
		std::for_each(objects_.begin(), objects_.end(), delete_ptr<HtmlObject> );
	}
}

//---------------------------------------------------------------------------------------------------------------------
std::ostream & operator << ( std::ostream & out, HtmlObject & obj)
{
	obj.render( out );
	return out;
}

//---------------------------------------------------------------------------------------------------------------------
std::ostream & operator << ( std::ostream & out, HtmlObject * obj)
{
	obj->render( out );
	return out;
}

//---------------------------------------------------------------------------------------------------------------------
std::ostream & operator << ( std::ostream & out, const std::list<HtmlObject *> & obj)
{
	std::copy(obj.begin(), obj.end(), std::ostream_iterator<HtmlObject *>(out,"\n") );
	return out;
}

//---------------------------------------------------------------------------------------------------------------------
HtmlHeader::HtmlHeader(const std::string & title, const std::string & author)
	: title_(title), author_(author)
{}

//---------------------------------------------------------------------------------------------------------------------
HtmlHeader::~HtmlHeader()
{}

//---------------------------------------------------------------------------------------------------------------------
RENDER(HtmlHeader,"head")

//---------------------------------------------------------------------------------------------------------------------
HtmlBody::HtmlBody() 
{}

//---------------------------------------------------------------------------------------------------------------------
HtmlBody::~HtmlBody()
{}

//---------------------------------------------------------------------------------------------------------------------
RENDER(HtmlBody,"body")

//---------------------------------------------------------------------------------------------------------------------
HtmlTag::HtmlTag(const std::string& tag)
	: tag_(tag)
{}

//---------------------------------------------------------------------------------------------------------------------
HtmlTag::~HtmlTag()
{}

//---------------------------------------------------------------------------------------------------------------------
void HtmlTag::render(std::ostream & out)
{
	if( text_.str().empty() && objects_.empty() ) {
		out << "<" << tag_ << attributes_ << "/>\n";
	} else {
		out << indent() << "<" << tag_ << attributes_ << ">"	
		    << text_.str() << "\n";			
		out << objects_;				
		out << unindent() << "</" << tag_ << ">";	
	}
}

//---------------------------------------------------------------------------------------------------------------------
HtmlPlot::HtmlPlot(TCanvas * canvas, bool remove, const std::string& name,bool pdf, bool macro, bool root)
	: canvas_(canvas), remove_(remove), name_(name), pdf_(pdf), macro_(macro), root_(root)
{
	if( name_.empty() && canvas_ !=0 ) { name_ = canvas_->GetName(); }
}

//---------------------------------------------------------------------------------------------------------------------
HtmlPlot::~HtmlPlot()
{}

//---------------------------------------------------------------------------------------------------------------------
void HtmlPlot::caption(const std::string & x) { caption_ =  x; } 


//---------------------------------------------------------------------------------------------------------------------
void HtmlPlot::render(std::ostream & out)
{
	if( canvas_ != 0 ) {
		canvas_->SaveAs((name_+".png").c_str());
		if( pdf_ ) {
			canvas_->SaveAs((name_+".pdf"  ).c_str());
		}
		if ( macro_ ) {
			canvas_->SaveAs((name_+".C"  ).c_str());
		}
		if ( root_ ) {
			canvas_->SaveAs((name_+".root"  ).c_str());
		}
	}
	if( remove_ ) { 
		delete canvas_; 
		canvas_ = 0;
	}
	out << indent() << "<div " << attributes_ << ">\n";
	if(!caption_.empty()) {
		out  << "<p>" << caption_ << /*"</p>" <<*/ "\n";
	}
	if( pdf_ ) {
		out << "<a href=\""  << name_ << ".pdf\"" << "> pdf" << "</a>";
	}
	if( macro_ ) {
		out << "<a href=\""  << name_ << ".C\"" << "> C" << "</a>";
	}
	if( root_ ) {
		out << "<a href=\""  << name_ << ".root\"" << "> root" << "</a>";
	}
	out << "<a href=\""  << name_ << ".png\"" << "> png" << "</a>";
	out << "<img src=\"" << name_ << ".png" << "\" />\n" ;
	out << "</p>" << "\n";

	out << unindent() << "</div>";
}

//---------------------------------------------------------------------------------------------------------------------
HtmlTable::HtmlTable()
{
	tbody_ = new HtmlTag("tbody");
	add(tbody_);
}

//---------------------------------------------------------------------------------------------------------------------
HtmlTable::~HtmlTable()
{}

//---------------------------------------------------------------------------------------------------------------------
HtmlTable::Row & HtmlTable::row()
{
	return tbody().add<Row>(new Row());
}

//---------------------------------------------------------------------------------------------------------------------
RENDER(HtmlTable,"table")

//---------------------------------------------------------------------------------------------------------------------
HtmlTable::Cell::Cell(HtmlObject *obj) 
{
	if( obj ) { add(obj); }
}

//---------------------------------------------------------------------------------------------------------------------
HtmlTable::Cell::~Cell()
{}

//---------------------------------------------------------------------------------------------------------------------
RENDER(HtmlTable::Cell,"td")

//---------------------------------------------------------------------------------------------------------------------
HtmlTable::Row::Row()
{}

//---------------------------------------------------------------------------------------------------------------------
HtmlTable::Row::~Row()
{}

//---------------------------------------------------------------------------------------------------------------------
HtmlTable::Cell & HtmlTable::Row::cell(HtmlObject * obj)
{
	return add<Cell>(new Cell(obj));
}

//---------------------------------------------------------------------------------------------------------------------
RENDER(HtmlTable::Row,"tr")

//---------------------------------------------------------------------------------------------------------------------
HtmlHelper::HtmlHelper(const std::string& dirname, const std::string &fname )
	: dirname_(dirname),
	  filename_(fname)
{
	add( &header_ );
	add( &body_ );
	navbar_ = &body().add(new HtmlTable()).row();
	own(false);
}

//---------------------------------------------------------------------------------------------------------------------
HtmlHelper::~HtmlHelper()
{}

//---------------------------------------------------------------------------------------------------------------------
void dumpIt(HtmlHelper * o)
{
	o->dump(false);
}


//---------------------------------------------------------------------------------------------------------------------
void HtmlHelper::dump(bool mkdir)
{
	if( mkdir ) {
		gSystem->mkdir(dirname_.c_str());
		gSystem->cd(dirname_.c_str());
	}
	std::ofstream fout(filename_.c_str(),std::ios_base::trunc);
	
	fout << *this << std::endl;
	fout.close();
	std::for_each(subpages_.begin(),subpages_.end(),dumpIt);
}

//---------------------------------------------------------------------------------------------------------------------
HtmlHelper & HtmlHelper::addPage(const std::string& name)
{
	navbar().cell().add<HtmlTag>(new HtmlTag("a")).set("href",name+".html") << name;
	HtmlHelper * h = new HtmlHelper(dirname_,name+".html");
	h->body().firstChild().add( &navbar() );
	subpages_.push_back(h);
	return **subpages_.rbegin();
}

//---------------------------------------------------------------------------------------------------------------------
RENDER(HtmlHelper,"html")

// Local Variables:
// mode: c++
// mode: sensitive
// c-basic-offset: 8
// End:

