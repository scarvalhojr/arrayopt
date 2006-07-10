#ifndef _SORT_H
#define _SORT_H


#define define_heap_operation(my_type) \
void upheap_##my_type(my_type a[], void *s, int n) \
{ \
  int i = n/2; \
  while (n > 1 && heap_compare(n, i)) { \
    heap_swap(i,n);\
    n = i; \
    i = n/2; \
  } \
} \
void heapfy_##my_type(my_type a[], void *s, int i, int n) \
{ \
  int j; \
  while (i <= n/2){ \
    j = 2*i; \
    if (j<n && heap_compare((j+1), j)) j++; \
    if (heap_compare(i,j)) break; \
    heap_swap(i, j); \
    i = j; \
  } \
} \
void heap_sort_##my_type(my_type a[], void *s, int n) \
{ \
  int i; \
  for (i=2; i<=n; i++) \
    upheap_##my_type (a, s, i); \
  for (i=n; i>1; i--){ \
    heap_swap(1, i); \
    heapfy_##my_type (a, s, 1, i-1); \
  } \
} \
my_type extract_first_##my_type(my_type a[], void *s, int n) \
{ \
  heap_swap(1, n); \
  heapfy_##my_type(a, s, 1, n-1); \
  return a[n]; \
} \
void heap_adjust_##my_type(my_type a[], void *s, int i, int n) \
{ \
  upheap_##my_type (a, s, i); \
  heapfy_##my_type (a, s, i, n); \
}
 

#define define_extended_heap_operation(my_type) \
void upheap_##my_type(my_type a[], void *s, int n) \
{ \
  int i = n/2; \
  while (n > 1 && heap_compare(n, i)) { \
    heap_swap(i,n);\
    n = i; \
    i = n/2; \
  } \
} \
void heapfy_##my_type(my_type a[], void *s, int i, int n) \
{ \
  int j; \
  while (i <= n/2){ \
    j = 2*i; \
    if (j<n && heap_compare((j+1), j)) j++; \
    if (heap_compare(i,j)) break; \
    heap_swap(i, j); \
    i = j; \
  } \
} \
void heap_sort_##my_type(my_type a[], void *s, int n) \
{ \
  int i; \
  for (i=2; i<=n; i++) \
    upheap_##my_type (a, s, i); \
  for (i=n; i>1; i--){ \
    heap_swap(1, i); \
    heapfy_##my_type (a, s, 1, i-1); \
  } \
} \
void heap_sortn_##my_type(my_type a[], void *s, int n, int k) \
{ \
   int i; \
   for (i=2; i<=n; i++) \
	 upheap_##my_type (a, s, i); \
   for (i=k; i>1; i--){ \
       heap_swap(1, i); \
       heapfy_##my_type (a, s, 1, i-1); \
    } \
} \
my_type extract_first_##my_type(my_type a[], void *s, int n) \
{ \
  heap_swap(1, n); \
  heapfy_##my_type(a, s, 1, n-1); \
  return a[n]; \
} \
void heap_adjust_##my_type(my_type a[], void *s, int i, int n) \
{ \
  upheap_##my_type (a, s, i); \
  heapfy_##my_type (a, s, i, n); \
}
 

#endif /* _SORT_H */
