clear; clc;
load("Mdata\QuadrotorDataTorqueDis12.mat");


% Export 3D array Dinput0 (3×30×30) to Excel with blank row separators
% Each Dinput0(:,:,k) is a 3×30 matrix, 30 blocks total
% Blank rows are inserted between blocks

%% File and data parameters
outputFilename = 'Dinput0_exp.xlsx';  % Output filename
sheetName = 'Sheet1';                     % Target sheet name
blankRows = 3;                            % Number of blank rows between blocks

%% Initialize combined data matrix
% Calculate total rows: 
% 3 rows per block × 30 blocks + 3 blank rows × 29 gaps
totalRows = size(Dinput0,1)*size(Dinput0,3) + blankRows*(size(Dinput0,3)-1);
combinedData = nan(totalRows, size(Dinput0,2));  % Preallocate with NaN

%% Fill data with blank row separators
currentRow = 1;  % Track position in combined matrix

for blockNum = 1:size(Dinput0, 3)
    % Insert current block (3×30)
    combinedData(currentRow:currentRow+size(Dinput0,1)-1, :) = Dinput0(:,:,blockNum);
    currentRow = currentRow + size(Dinput0,1);
    
    % Insert blank rows (except after last block)
    if blockNum < size(Dinput0, 3)
        currentRow = currentRow + blankRows;
    end
end

%% Write to Excel
writematrix(combinedData, outputFilename, 'Sheet', sheetName);

%% Verification
% Check if file was created
if isfile(outputFilename)
    fprintf('Success! Data exported to: %s\n', outputFilename); 
    fprintf('Total data rows: %d (including %d blank rows)\n', ...
            size(combinedData,1), blankRows*(size(Dinput0,3)-1));
else
    error('File creation failed. Check write permissions.');
end

































