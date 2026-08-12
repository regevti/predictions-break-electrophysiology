function [preMask, postMask, preDuration, postDuration] = prediction_break_epoch_masks(t, options)
%PREDICTION_BREAK_EPOCH_MASKS Shared half-open pre/post event windows.
%   PRE is [preWindow(1), preWindow(2)); POST is
%   [postWindow(1), postWindow(2)). Half-open windows assign t=0 to POST
%   and prevent the event boundary from being counted twice.

    arguments
        t double
        options.preWindow (1,2) double = [-1 0]
        options.postWindow (1,2) double = [0 1]
    end

    validateattributes(options.preWindow, {'double'}, {'increasing'});
    validateattributes(options.postWindow, {'double'}, {'increasing'});

    % Colon-generated bin starts can represent zero as a tiny negative
    % number (for example -1.11e-16). Shift comparisons by a tolerance far
    % below the sampling interval so that the nominal zero bin belongs to
    % POST, not PRE.
    tolerance = 100 * eps(max(1,max(abs([t(:); options.preWindow(:); options.postWindow(:)]))));
    preMask = t >= options.preWindow(1)-tolerance & t < options.preWindow(2)-tolerance;
    postMask = t >= options.postWindow(1)-tolerance & t < options.postWindow(2)-tolerance;
    preDuration = diff(options.preWindow);
    postDuration = diff(options.postWindow);

    if ~any(preMask) || ~any(postMask)
        error('predictionBreak:EmptyEpoch', ...
              'The time vector does not contain samples in both requested epochs.');
    end
    if any(preMask & postMask)
        error('predictionBreak:OverlappingEpochs', ...
              'Pre and post masks overlap; use non-overlapping half-open windows.');
    end
end
