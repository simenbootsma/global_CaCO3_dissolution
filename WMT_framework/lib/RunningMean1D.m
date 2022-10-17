function [smooth] = RunningMean1D(vector,elements)
% A function to smooth a 2D matrix,.
% In particular tracer gradients.

%% INPUT
% Matrix        = 1xN vector that is to be smoothed.
% bins_columns  = Number of values to average in the N-direction (in a column).
% repeats       = The number of time we repeat this operation.

%% OUTPUT
% A Matrix of NxM length smoothes.
% Each (n,m) point is now an average ofOutput

% Reshape matrix
[M,N]       = size(vector);
if M>1
    vector = vector';
    [~,N]  = size(vector);
end
NaNmask     = vector-vector;
smooth      = NaN(1,N);

    for n = 1:N % Number of elements to average.
            nb = n-elements; % Start index column.
            ne = n+elements; % End index column
            % If statements for if we hit the start or end of the domain:
            if nb <= 0  % If we hit the edge of the domain.
                nb=1;
            end
            if ne >= N % If we hit the end of the domain.
                ne = N;
            end

            % Select values to be averaged:
            val             = vector(nb:ne);
            % Remove the NaNs
            val(isnan(val)) = [];
            % Average them:
            if isempty(val) % set to NaN, if no values.
               smooth(n)=NaN;
            else % Average otherwise:
               smooth(n) = mean(val(:));
            end
            clear val ne nb
    end

% Final matrix we want the NaNs back:
smooth = smooth + NaNmask;