#ifndef __LBM_LOG_H__
#define __LBM_LOG_H__

#include <stdio.h>
#include <time.h>


#ifdef NDEBUG
#define DEBUG(M, ...)
#else
#define DEBUG(M, ...) fprintf(stderr, "[%s  %s] [DEBUG] (%d): " M "", __DATE__, __TIME__, __LINE__, ##__VA_ARGS__)
#endif

#ifdef NMLOG

#define MLOG(M, ...)

#else

#define MLOG(M, ...) {if(myrank == 0) { \
	                  time_t t; \
					  struct tm *tmif; \
	                  t = time(NULL); \
	                  tmif = localtime(&t); \
	                  fprintf(stderr, "[%d-%2d-%2d] [%2d:%2d:%2d] [MLOG] " M "",  \
							  tmif->tm_year + 1900, \
							  tmif->tm_mon + 1, \
							  tmif->tm_mday, \
							  tmif->tm_hour, \
							  tmif->tm_min, \
							  tmif->tm_sec, ##__VA_ARGS__); }} 
#endif

#define TLOG(M, ...)  {if(myrank == 0) { \
	                  fprintf(stderr, "[%s  %s] [TLOG ] (%d): " M "", __DATE__, __TIME__, __LINE__, ##__VA_ARGS__) }}

#define OLOG(RANK, M, ...)  \
	                  fprintf(stderr, "[%s  %s] [OLOG ] (RANK%d  %d): " M "", __DATE__, __TIME__, RANK, __LINE__, ##__VA_ARGS__)


#endif
