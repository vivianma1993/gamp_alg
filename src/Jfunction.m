function [J] = Jfunction(alpha)
% JFunction is an approximated function for EXIT curves of AWGN
% channel model, it takes alpha vector and returns a corresponding
% JFunction vector.
% Values are arbitary that fullfill required function
alpha_tp = 1.6363;
a_J1 = -0.0421061;
b_J1 = 0.209252;
c_J1 = -0.00640081;
a_J2 = 0.00181491;
b_J2 = -0.142675;
c_J2 = -0.0822054;
d_J2 = 0.0549608;
J_temp = zeros(1,length(alpha));

J_temp(alpha >= 0 & alpha <= alpha_tp) = (a_J1*(alpha(alpha >= 0 & alpha <= alpha_tp).^3))+(b_J1*(alpha(alpha >= 0 & alpha <= alpha_tp).^2))+(c_J1*(alpha(alpha >= 0 & alpha <= alpha_tp)));
J_temp(alpha > alpha_tp & alpha < 10) = 1-exp((a_J2*(alpha(alpha > alpha_tp & alpha < 10).^3))+(b_J2*(alpha(alpha > alpha_tp & alpha < 10).^2))+(c_J2*(alpha(alpha > alpha_tp & alpha < 10)))+(d_J2));
J_temp(alpha >= 10) = 1;

J = J_temp;
end
