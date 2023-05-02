function t_w = warp_t(t,r,c)
%%% warp_t.m
%%% Julien Bonnel, Woods Hole Oceanographic Institution
%%% March 2019

%%% warping subroutine

t_w=sqrt(t.^2+r^2/c^2);
