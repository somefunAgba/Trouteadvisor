function new_path = move_sets(NumCity,path)
% set new tour path to path
new_path = path;
%% MOVE: RouletteWheelSelection:
p_invert = 0.5; % inversion probability
p_swap = 0.2; % switch probability
p_slide = 1-p_invert-p_swap; % translation probability
rnum = rand;
csum=cumsum([p_invert p_swap p_slide]);
move=find(rnum<=csum,1,'first');

switch move
    case 1
        % FLIP: perform an INVERSION: randomly remove a section of the tour
        % specified by two endpoints and replace them in the opposite order.
        % e.g 6789 becomes 9876
        index = randi([2,NumCity],2,1);
        inversion_index = (min(index):max(index));
        % invert the new tour path
        new_path(inversion_index) = fliplr(path(inversion_index));
    case 2
        % SWAP: perform a SWITCH: randomly select two endpoints,
        % cities and switch them in a tour e.g 3 11 becomes 11 3
        index = randsample([2,NumCity],2);
        swap_index = [min(index) max(index)];
        new_path(swap_index) = path([max(index) min(index)]);
    case 3
        % SLIDE: perform a TRANSLATION: remove a section of the tour specified by two
        % endpoints and replace it between two randomly selected consecutive cities
        % e.g: remove 8 7 and slide between random 4 5
        index = randsample([2,NumCity],2,1);
        slide_index = (min(index):max(index));
        new_path(slide_index) = path([min(index)+1:max(index) min(index)]);
end
end