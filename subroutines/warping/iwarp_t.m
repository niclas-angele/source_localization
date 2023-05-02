function [ t ] = iwarp_t(t_w,r,c)
%%% iwarp_t.m
%%% Julien Bonnel, Woods Hole Oceanographic Institution
%%% March 2019

%%% warping subroutine

t=sqrt(t_w.^2-(r/c)^2);