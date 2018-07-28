#ifndef THREADS_H__
#define THREADS_H__
#include "cmdlib.h" //--vluzacn

#include <pointerTo.h>
#include <lightMutex.h>
#include <genericThread.h>
#include <threadPriority.h>
#include <pvector.h>

#if _MSC_VER >= 1000
#pragma once
#endif

#define	MAX_THREADS	64

typedef void q_threadfunction( int );

#ifdef SYSTEM_WIN32
#define DEFAULT_NUMTHREADS -1
#endif
#ifdef SYSTEM_POSIX
#define DEFAULT_NUMTHREADS 1
#endif

class BSPThread : public Thread
{
public:
        BSPThread();
        void set_function( q_threadfunction *func );
        void set_value( int val );
        volatile bool is_finished() const;

protected:
        virtual void thread_main();

private:
        q_threadfunction * _func;
        int _val;
        volatile bool _finished;
};

#define DEFAULT_THREAD_PRIORITY TP_normal

extern int      g_numthreads;
extern ThreadPriority g_threadpriority;
extern LightMutex g_global_lock;
extern pvector<PT( BSPThread )> g_threadhandles;

extern int	GetCurrentThreadNumber();

extern void     ThreadSetPriority( ThreadPriority type );
extern void     ThreadSetDefault();
extern int      GetThreadWork();
extern void     ThreadLock();
extern void     ThreadUnlock();

extern void     RunThreadsOnIndividual( int workcnt, bool showpacifier, q_threadfunction );
extern void     RunThreadsOn( int workcnt, bool showpacifier, q_threadfunction );

#ifdef ZHLT_NETVIS
extern void     threads_InitCrit();
extern void     threads_UninitCrit();
#endif

#define NamedRunThreadsOn(n,p,f) { Log("%s\n", #f); RunThreadsOn(n,p,f); }
#define NamedRunThreadsOnIndividual(n,p,f) { Log("%s\n", #f); RunThreadsOnIndividual(n,p,f); }

#endif //**/ THREADS_H__
