#ifndef HtmlHelper_h
#define HtmlHelper_h

#include <string>
#include <sstream>

#include <list>
#include <map>


/**
 *
 *
 */
class HtmlObject
{
public:
	HtmlObject();
	virtual ~HtmlObject();
	
	virtual void render(std::ostream &) = 0;
	
	HtmlObject & set(const std::string&, const std::string&);
	
	template<class T> T & add(T * obj) { 
		objects_.push_back((HtmlObject*)obj);  
		return *obj; 
	}
	
	HtmlObject & addobj(HtmlObject * obj) {
             objects_.push_back((HtmlObject*)obj);
             return *obj;
        }
	
	template<class T> HtmlObject & operator << (const T & t) {
		text_ << t;
		return *this;
	}

	HtmlObject & txt(const std::string & st) {
		*this << st;
		return *this;
	}

	void own(bool x) { own_ = x; }
	bool own()       { return own_; }

	HtmlObject & firstChild() { return objects_.empty() ? *this : *(*objects_.begin());  } 
	
protected:
	std::string indent();
	std::string unindent();

	static int indent_;
	
	bool own_;
	std::stringstream text_;
	std::list<HtmlObject *> objects_;
	std::list<std::pair<std::string,std::string> > attributes_;
};

std::ostream & operator << ( std::ostream &, HtmlObject & );
std::ostream & operator << ( std::ostream &, HtmlObject * );

std::ostream & operator << ( std::ostream &, const std::list<HtmlObject *> & );

/**
 *
 *
 */
class HtmlHeader : public HtmlObject 
{
public:
	HtmlHeader(const std::string & title="", const std::string & author="");
	virtual ~HtmlHeader();
	
	void render(std::ostream &);
private:
	std::string title_;
	std::string author_;
};

/**
 *
 *
 */
class HtmlBody : public HtmlObject 
{
public:
	HtmlBody();
	virtual ~HtmlBody();
	
	void render(std::ostream &);
};

/**
 *
 *
 */
class HtmlTag : public HtmlObject
{
public:
        HtmlTag(const std::string& tag);
        virtual ~HtmlTag();

        void render(std::ostream &);
	
private:
	std::string tag_;
	
};

/**
 *
 *
 */
class TCanvas;
class HtmlPlot : public HtmlObject {
public:
	HtmlPlot(TCanvas *, bool remove=false, const std::string& name = "", bool pdf=true, bool macro=false, bool root=false);
	virtual ~HtmlPlot();
	
	void caption(const std::string &);
	
	void render(std::ostream &);

private:
	TCanvas * canvas_;
	bool remove_;
	std::string name_; 
	std::string caption_;	
	bool pdf_, macro_, root_;
	
};

/**
 *
 *
 */
class HtmlTable : public HtmlObject 
{
public:
	HtmlTable();
	virtual ~HtmlTable();
	
	void render(std::ostream &);
	
	class Cell  : public HtmlObject {
	public:
		Cell(HtmlObject * obj);
		virtual ~Cell();
		
		void render(std::ostream &);
	};
	
	class Row  : public HtmlObject {
	public:
		Row();
		virtual ~Row();
		
		Cell & cell(HtmlObject * obj=0); 
		void render(std::ostream &);
	};
	
	Row & row();

	HtmlTag & tbody() { return * tbody_; } 
	
private:
	HtmlTag * tbody_;
};

/**
 *
 *
 */
class HtmlHelper : public HtmlObject 
{
public:
	HtmlHelper(const std::string& dirname, const std::string& fname="index.html");
	virtual ~HtmlHelper();
	
	void dump(bool mkdir=true);
	
	HtmlHeader & header() { return header_; }
	HtmlBody   & body()   { return body_  ; }
	HtmlTable::Row  & navbar() { return *navbar_; }
	
	HtmlHelper & addPage(const std::string& name);
	
	void render(std::ostream &);
	
private:
	std::string dirname_; 
	std::string filename_; 
	HtmlHeader header_;
	HtmlBody   body_;
	HtmlTable::Row  * navbar_;
	
	std::list<HtmlHelper*> subpages_;
};


template HtmlHeader     & HtmlObject::add<HtmlHeader     >(HtmlHeader*);
template HtmlBody       & HtmlObject::add<HtmlBody       >(HtmlBody  *);
template HtmlPlot       & HtmlObject::add<HtmlPlot       >(HtmlPlot  *);
template HtmlTable      & HtmlObject::add<HtmlTable      >(HtmlTable *);
template HtmlTable::Row & HtmlObject::add<HtmlTable::Row >(HtmlTable::Row *);
template HtmlTable::Cell& HtmlObject::add<HtmlTable::Cell>(HtmlTable::Cell *);
template HtmlHelper     & HtmlObject::add<HtmlHelper     >(HtmlHelper*);
template HtmlObject     & HtmlObject::operator<< (const std::string  & );

#endif

// Local Variables:
// mode: c++
// mode: sensitive
// c-basic-offset: 8
// End:

