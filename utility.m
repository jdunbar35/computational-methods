%% Q6: Household utility funciton
% Jack Dunbar
% October 31, 2024

function utility = utility(c, g, l, gamma)
    utility = log(c) + gamma * log(g) - l^2 / 2;
end