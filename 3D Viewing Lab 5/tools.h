#pragma once
#pragma once
#include<iostream>
#include<math.h>

#include"structures.h"

void bubbleSort(vertex(&arr)[3], int n)
{
    int i, j;
    for (i = 0; i < n - 1; i++)

        // Last i elements are already in place 
        for (j = 0; j < n - i - 1; j++)
            if (arr[j].y > arr[j + 1].y)
                swap(arr[j], arr[j + 1]);

}

int partition(vertex(&arr)[3], int low, int high)
{
    int pivot = arr[high].y;    // pivot
    int i = (low - 1);  // Index of smaller element

    for (int j = low; j <= high - 1; j++)
    {
        // If current element is smaller than or
        // equal to pivot
        if (arr[j].y <= pivot)
        {
            i++;    // increment index of smaller element
            swap(arr[i], arr[j]);
        }
    }
    swap(arr[i + 1], arr[high]);
    return (i + 1);
}

void quickSort(vertex(&arr)[3], int low, int high)
{
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now
           at right place */
        int pi = partition(arr, low, high);

        // Separately sort elements before
        // partition and after partition
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

void swap(vertex* xp, vertex* yp)
{
    vertex temp = *xp;
    *xp = *yp;
    *yp = temp;
}

//scanline algorithm

//void RasterizeTriangle(vertex v1, vertex v2, vertex v3) {
//    vertex arr[3] = { v1,v2,v3 };
//   // v1.y < v2.y ? swap(v2, v1) : v1.y < v3.y ? swap(v3, v1) : v2.y < v3.y ? swap(v3, v2) : ;
//    arr[0].y < arr[1].y ? swap(arr[1], arr[0]) : arr[0].y < arr[2].y ? swap(arr[2], arr[0]) : arr[1].y < arr[2].y ? swap(arr[2], arr[1]) :swap(arr[1],arr[1]) ;
//    if (v1.y = v2.y)return;
//
//
//
//
//
//}