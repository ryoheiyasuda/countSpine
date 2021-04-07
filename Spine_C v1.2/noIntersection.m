function [ bool ] = noIntersection( highlightImageA, highlightImageB )
%NOINTERSECTION Summary of this function goes here
%   Detailed explanation goes here
combination = highlightImageA + highlightImageB;
bool = ~any(combination(:) == 2);