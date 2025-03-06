#include "dd/DDDebug.hpp"
#include <cstdarg>

char msgBuffer[MSG_BUFFER_LENGHTH];

static void putInt(int n) {
  if(n > 9) {
    putInt(n/10);
  }
  putchar(n % 10 + '0');
}

int debug_info_printf(const char *format, ...) {
  char buffer[MSG_BUFFER_LENGHTH];
  char *cursor = buffer;

  if(strlen(format) == 0) {
    return -1;
  }
  va_list arg;
  va_start(arg, format);

  const char *start = format;

  while(*start!='\0') {
    if(*start == '%') {
      start += 1;
      switch (*start) {
        case 'd':
          sprintf(cursor, "%d", va_arg(arg, int));
          cursor = buffer + strlen(buffer);
          break;
        case 's':
          char *ch = va_arg(arg, char*);
          char char_buf[128];
          char *c = char_buf;
          while(*ch) {
            strcpy(c, ch);
            c++;
            ch++;
          }
          sprintf(cursor, "%s", char_buf);
          cursor = buffer + strlen(buffer);
          break;
      }
    } else {
      sprintf(cursor, "%c", *start);
      cursor = buffer + strlen(buffer);
    }
    ++start;
  }
  va_end(arg);
  DEBUG_INFO(buffer);

  return 0;
}

int debug_error_printf(const char *format, ...) {
  char buffer[MSG_BUFFER_LENGHTH];
  char *cursor = buffer;

  if(strlen(format) == 0) {
    return -1;
  }
  va_list arg;
  va_start(arg, format);

  const char *start = format;

  while(*start!='\0') {
    if(*start == '%') {
      start += 1;
      switch (*start) {
        case 'd':
          sprintf(cursor, "%d", va_arg(arg, int));
          cursor = buffer + strlen(buffer);
          break;
        case 's':
          char *ch = va_arg(arg, char*);
          char char_buf[128];
          char *c = char_buf;
          while(*ch) {
            strcpy(c, ch);
            c++;
            ch++;
          }
          sprintf(cursor, "%s", char_buf);
          cursor = buffer + strlen(buffer);
          break;
      }
    } else {
      sprintf(cursor, "%c", *start);
      cursor = buffer + strlen(buffer);
    }
    ++start;
  }
  va_end(arg);
  DEBUG_ERROR(buffer);

  return 0;
}
