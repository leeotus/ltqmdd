#include <stdio.h>

char *impl_strcpy(char *dst, const char *src) {
  char *res = dst;
  while((*dst++ = *src++) != '\0');
  return res;
}

int impl_strlen(const char *str) {
  int len = 0;
  while((*str++) != '\0') {
    ++len;
  }
  return len;
}

char *impl_strcat(char *dst, const char *src) {
  char *it = dst;
  while((*it) != '\0') {
    it++;
  }
  while((*it++ = *src++) != '\0');
  return dst;
}

int impl_strcmp(const char *str1, const char *str2) {
  while(*str1 && *str2 && (*str1 == *str2)) {
    str1++;
    str2++;
  }
  return *str1 - *str2;
}

int main(int argc, char **argv) {
  char buf[20];
  char *hello = "hello ";
  char *world = "world";
  char *b = impl_strcat(buf, hello);
  b = impl_strcat(buf, world);
  printf("%s\r\n", b);
  printf("%d\r\n", impl_strlen(b));
  return 0;
}
