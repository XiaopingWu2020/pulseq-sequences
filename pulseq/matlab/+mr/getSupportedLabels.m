function supported_labels = getSupportedLabels()
% auxilary function

supported_labels={'SLC','SEG','REP','AVG','SET','ECO','PHS','LIN','PAR',... % labels/counters
    'NAV','REV','SMS','PMC','DUM','RTF',... % flags
    'GT','RO'}; % flags not related to adc that are used to identify gradient blocks that are skope relevant.

end
