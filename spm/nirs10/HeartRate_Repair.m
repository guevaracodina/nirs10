function [OutRate GapFlag TotalAdded] = HeartRate_Repair(InRate, Precision, Gap_def)
%**************************************************************************
%Function: this function detects the gaps between heart rate records and 
%          fill the gaps automatically.
%          InRate: Original Heart Rate sets
%          Precision: Int value, wide-range average gap detection or 
%                     narrow-range average gap detection
%          Gap_def: Double value, indicates the definition of gaps. For one
%                   number in the set, how many times of the average value 
%                   could it be recognised as a gap value. However, the 
%                   number of gaps to be divided will use the round()
%                   calculation.
%          OutRate: Repaired Heart Rate sets
%          GapFlag: Int value, indicates whether there is a gap in the 
%                   original sets (0-No, 1-Yes)
%          TotalAdded: Int value, indicates how many values in total has
%                      been added to the repaired heart rate sets
%Editor: Ke Peng, Laboratoire d'Imagerie Optique et Moleculaire
%Data: March 08, 2012
%**************************************************************************

InRate_order = InRate;

check_queue = zeros(1,Precision);
GapFlag = 0;             %Flag indicating whether there is a gap in between
TotalAdded = 0;         %indicating how many data has been added in between

if length(InRate) <= Precision
       %If total length of HeartRate record is less than precision demanded
    TotalRate = 0;
    for i_check = 1 : length(InRate)
        TotalRate = TotalRate + InRate(i_check);
    end
        
    for i_check = 1 : length(InRate)
        
        i_TotalRate = TotalRate - InRate(i_check);
        i_MeanRate = i_TotalRate/(length(InRate)-1);
        
        i_generate = i_check;
        i_generate = i_generate + TotalAdded;
        
        if InRate(i_check) < i_MeanRate*Gap_def
            continue;
        else
            GapFlag = 1;
            Gapnum = round(InRate(i_check)/i_MeanRate);
            Rate_add = InRate(i_check) / Gapnum;
            InRate_order(i_generate) = Rate_add;
            
            for i_add = 1 : Gapnum-1
                InRate_order = [InRate_order(1:i_generate) Rate_add ...
                    InRate_order(i_generate+1 : length(InRate_order))];
                TotalAdded = TotalAdded + 1;
            end
            
        end
    end
    
else

    for i_check = 1 : length(InRate)
        i_generate = i_check;
        i_generate = i_generate + TotalAdded;    
        if i_generate == 1                          %if the first HeartRate
            for i_queue = 1 : Precision
                check_queue(i_queue) = InRate_order(i_queue + 1);
            end
        elseif i_generate <= Precision         %if the first five HeartRate
            for i_queue = 1 : i_check-1
                check_queue(i_queue) = InRate_order(i_queue);
            end
            for j_queue = i_generate+1 : i_generate+Precision-i_queue
                check_queue(j_queue-1) = InRate_order(j_queue);
            end
            clear j_queue
        else                                            %if other HeartRate
            for i_queue = 1 : Precision
                check_queue(i_queue) = InRate_order(i_generate - i_queue);
            end
        end

        clear i_queue

        TotalRate = 0;
        for i_mean = 1:Precision
            TotalRate = TotalRate + check_queue(i_mean);
        end
        MeanRate = TotalRate / Precision;
        clear TotalRate

        if InRate(i_check) < MeanRate*Gap_def
            continue;
        else
            GapFlag = 1;                                         %outputted
            Gapnum = round(InRate(i_check) / MeanRate);
            Rate_add = InRate(i_check) / Gapnum;

            InRate_order(i_generate) = Rate_add;
            for i_add = 1 : Gapnum-1
                InRate_order = [InRate_order(1:i_generate) Rate_add ...
                    InRate_order(i_generate+1 : length(InRate_order))];
                TotalAdded = TotalAdded+1;
            end
        end
    end
end

OutRate = InRate_order;

