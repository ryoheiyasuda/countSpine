function [] = spineExcel()
cs = get(gcf, 'UserData');
spineNum = length(cs.data.spineLength);
fileName = uiputfile('*.xls','Save Data As...');
numbers(1) = {'Number'};
for i = 2:spineNum + 1;
    numbers(i) = {i - 1};
end
numbersCells = transpose(numbers);
spineLength = cs.data.spineLength;
spineLengthCells(1) = {'Spine Length(um)'};
for i = 1:spineNum
    spineLengthCells(i+1) = {spineLength(i)};
end
spineType = cs.data.spineType;
spineTypeCells(1) = {'Spine Type'};
for i = 1:spineNum
    spineTypeCells(i+1) = spineType(i);
end
dendLengthCells(1) = {'Dendrite Length(um)'};
dendLengthCells(2) = {cs.data.dendLength};
for i = 3:spineNum + 1
    dendLengthCells(i) = {NaN};
end
spineDensityCells(1) = {'Spine Density(spines/um)'};
spineDensityCells(2) = {cs.data.spineDensity};
for i = 3:spineNum + 1
    spineDensityCells(i) = {NaN};
end
dendLengthCells = transpose(dendLengthCells);
spineDensityCells = transpose(spineDensityCells);
spineLengthCells = transpose(spineLengthCells);
spineTypeCells = transpose(spineTypeCells);
excelMatrix = [numbersCells,spineLengthCells,spineTypeCells,dendLengthCells,spineDensityCells];
xlswrite(fileName, excelMatrix);
end