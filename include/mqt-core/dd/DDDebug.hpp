#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef enum __debug_type {
  WARNING,
  INFO,
  ERROR,
} DebugType;

typedef enum __colors {
  NONE,
  BLACK,
  L_BLACK,
  RED,
  L_RED,
  GREEN,
  L_GREEN,
  BROWN,
  YELLOW,
  BLUE,
  L_BLUE
} COLOR;

#define COLOR_NONE    "\e[0m"
#define COLOR_BLACK   "\e[0;30m"
#define COLOR_L_BLACK "\e[1;30m"
#define COLOR_RED     "\e[0;31m"
#define COLOR_L_RED   "\e[1;31m"
#define COLOR_GREEN   "\e[0;32m"
#define COLOR_L_GREEN "\e[1;32m"
#define COLOR_BROWN   "\e[0;33m"
#define COLOR_YELLOW  "\e[1;33m"
#define COLOR_BLUE    "\e[0;34m"
#define COLOR_L_BLUE  "\e[1;34m"
#define COLOR_CLOSE_TAG "\033[0m"

#define MSG_BUFFER_LENGHT 1024
extern char msgBuffer[MSG_BUFFER_LENGHT];
#define __formatted_msg(buffer, msg)                                           \
  do {                                                                         \
    bzero(&(buffer), MSG_BUFFER_LENGHT);                                       \
    sprintf(buffer, "In file %s, line %d, `%s`", (char*)__FILE__,              \
            (int)__LINE__, (char*)(msg));                                      \
  } while (0)

#define ENABLE_DYN_DEBUG
#if (defined(ENABLE_DYN_DEBUG))

#define DEBUG_WARNING(msg)                                                     \
  do {                                                                         \
    __formatted_msg(msgBuffer, msg);                                           \
    __debug_printf(msgBuffer, WARNING, GREEN);                                 \
  } while (0)

#define DEBUG_ERROR(msg)                                                       \
  do {                                                                         \
    __formatted_msg(msgBuffer, msg);                                           \
    __debug_printf(msgBuffer, ERROR, RED);                                      \
  } while (0)

#define DEBUG_INFO(msg)                                                        \
  do {                                                                         \
    __formatted_msg(msgBuffer, msg);                                           \
    __debug_printf(msgBuffer, INFO, BLUE);                                      \
  } while (0)

#else

#define DEBUG_WARNING(msg)
#define DEBUG_INFO(msg);
#define DEBUG_ERROR_COND(cond, msg)
#define DEBUG_ERROR(msg)

#endif

static void __debug_printf(void *msg, DebugType type, COLOR c)
{
  const char *color;
  switch (c) {
  case BLUE:
    color = (const char*)COLOR_BLUE;
    break;
  case RED:
    color = (const char*)COLOR_RED;
    break;
  case GREEN:
    color = (const char*)COLOR_GREEN;
    break;
  default:
    break;
  }
  switch (type)
  {
  case WARNING:
    printf("%s[ WARNING ]: %s %s\r\n", color, (char*)msg,
           (char*)COLOR_CLOSE_TAG);
    break;
  case INFO:
    printf("%s[ INFOMATION ]: %s %s\r\n", color, (char*)msg,
           (char*)COLOR_CLOSE_TAG);
    break;
  case ERROR:
    printf("%s[ ERROR ]: %s %s\r\n", color, (char*)msg,
           (char*)COLOR_CLOSE_TAG);
    break;
  default:
    break;
  }
}

int debug_info_printf(const char *format, ...);

int debug_error_printf(const char *format, ...);
