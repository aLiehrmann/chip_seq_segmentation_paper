#include <Rcpp.h>
#include<bits/stdc++.h>
using namespace Rcpp;

double kadane(double* arr, int* start, int* finish, int n)
{
    // initialize sum, maxSum and
    double sum = 0, maxSum = DBL_MIN;
    int i;

    // Just some initial value to check
    // for all negative values case
    *finish = -1;

    // local variable
    int local_start = 0;

    for (i = 0; i < n; ++i)
    {
        sum += arr[i];
        if (sum < 0)
        {
            sum = 0;
            local_start = i + 1;
        }
        else if (sum > maxSum)
        {
            maxSum = sum;
            *start = local_start;
            *finish = i;
        }
    }

    // There is at-least one
    // non-negative number
    if (*finish != -1)
        return maxSum;

    // Special Case: When all numbers
    // in arr[] are negative
    maxSum = arr[0];
    *start = *finish = 0;

    // Find the maximum element in array
    for (i = 1; i < n; i++)
    {
        if (arr[i] > maxSum)
        {
            maxSum = arr[i];
            *start = *finish = i;
        }
    }
    return maxSum;
}

// [[Rcpp::export]]
List max_sum_rectangle(NumericMatrix M){

    int finalLeft, finalRight, finalTop, finalBottom;
    double maxSum = DBL_MIN;
    int left, right, i, start, finish;
    double temp[M.nrow()], sum;
    for (left = 0; left < M.ncol(); ++left)
    {
        memset(temp, 0, sizeof(temp));
        for (right = left; right < M.ncol(); ++right)
        {
            for (i = 0; i < M.nrow(); ++i)
                temp[i] += M(i,right);
            sum = kadane(temp, &start, &finish, M.nrow());
            if (sum > maxSum)
            {
                maxSum = sum;
                finalLeft = left;
                finalRight = right;
                finalTop = start;
                finalBottom = finish;
            }
        }
    }
    List res = List::create(
        _["index.min.log.lambda"] = finalTop+1,
        _["index.max.log.lambda"] = finalBottom+1,
        _["index.min.log.phi"] = finalLeft+1,
        _["index.max.log.phi"] = finalRight+1,
        _["max.sum"] = maxSum
    );
    return res;
}