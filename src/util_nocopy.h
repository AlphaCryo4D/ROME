#pragma once

//===============================================================================================================================
// By using this as a member or as a base class, the code will get an error message
// if an attempt is made to copy the object.  
//
// Use when copying doesn't make sense, because its use causes obscure bugs to become compile-time errors.
//
class NoCopy {
	int _uid;
	NoCopy(NoCopy const & rhs) = delete;				// Private and unimplemented so can't be used
	NoCopy& operator=(NoCopy const & rhs) = delete;
public:
	__declspec(noinline) NoCopy();						// noinline avoids a IC++ 2018 beta test bug
    __declspec(noinline) void changeUid();				// noinline avoids a IC++ 2018 beta test bug
	int uid() const { return _uid; }
	bool interesting() const;
};
