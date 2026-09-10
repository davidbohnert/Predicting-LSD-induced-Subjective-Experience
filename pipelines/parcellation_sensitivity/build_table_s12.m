function [internalTable, crossDrugTable] = build_table_s12(resultRoot)
%BUILD_TABLE_S12 Compatibility alias for the renumbered Table S13 builder.
% New calls should use BUILD_TABLE_S13. Outputs use table_S13_*.csv;
% previously saved table_S12_*.csv files are left unchanged.
[internalTable, crossDrugTable] = build_table_s13(resultRoot);
end
