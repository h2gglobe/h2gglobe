#ifndef _Sorters_h_
#define _Sorters_h_

// -----------------------------------------------------------------------------------------------------
template<class T, class U, class C=std::less<U> > class ClonesSorter
{
public:
	typedef U (T::*method_type) () const;
	
	ClonesSorter(TClonesArray *a, method_type m, C c=C()) : array_(a), method_(m), compare_(c) {};
	
	bool operator () (int ii, int jj) {
		T * ei = dynamic_cast<T*> ((*array_)[ii]);
		T * ej = dynamic_cast<T*> ((*array_)[jj]);
		return compare_( (*ei.*method_)() ,  (*ej.*method_)() );
	};
	
private:
	TClonesArray * array_;
	method_type method_;
	C compare_;
};

// -----------------------------------------------------------------------------------------------------
template<class T, class C=std::less<T> > class SimpleSorter
{
public:
	SimpleSorter(T *a, C c=C()) : array_(a), compare_(c) {};
	
	bool operator () (int ii, int jj) {
		return compare_( array_[ii], array_[ii] );
	};
	
private:
	T * array_;
	C compare_;
};

#endif
 
