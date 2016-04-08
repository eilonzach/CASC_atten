function [IP,note] = which_OBS(sta)
% [IP,note] = which_OBS(sta)
% IP is 'WHOI', 'LDEO', or 'SIO';
% note is [] or 'TRM'

IP = [];
note = [];

whoi = {'G03A','G30A','J06A',...
        'J23A','J28A','J29A',...
        'J30A','J31A','J37A','J38A','J39A',...
        'J45A','J46A','J47A','J48A',...
        'J52A','J53A','J54A','J55A',...
        'J61A','J63A','J67A','J68A',...
        'FS05B','FS06B','FS09B',...
        'G03B','G04B','G05B',...
        'G11B','G13B','G19B',...
        'G20B','G21B','G22B','G29B',...
        'G30B','G35B','G36B2',...
        'J06B','J11','J19',...
        'J23','J27','J28',...
        'J48','J63',...
        'J21C','J23C','J28C','J29C',...
        'J30C','J31C','J32C','J35C','J36C','J37C','J38C','J39C',...
        'J43C','J44C','J45C','J46C','J47C','J48C',...
        'J54C','J55C','J63C','J67C','J69C',...
        };
        
ldeo = {'FN01A','FN05A','FN07A','FN08A',...
        'FN12A','FN14A','FN16A','FN18A',...
        'J26A','J34A','J41A','J42A','J49A',...
        'J50A','J51A','J58A','J59A',...
        'M03A','M06A',...
        'FS02B','FS03B','FS04B','FS07B','FS08B',...
        'FS10B','FS11B','FS12B','FS13B','FS15B','FS16B','FS17B','FS18B','FS19B',...
        'FS20B'...
        'G09B','G17B','G18B','G25B','G26B',...
        'G33B','G34B',...
        'J17B','J25B','J33B','J34B',...
        'M09B','M10B','M13B','M18B',...
        'FN01C','FN02C','FN03C','FN04C','FN05C','FN06C','FN07C','FN08C','FN09C',...
        'FN10C','FN11C','FN12C','FN13C','FN14C','FN16C','FN17C','FN18C','FN19C',...
        'J26C','J34C','J41C','J42C','J49C',...
        'J50C','J51C','J57C','J58C','J59C','M06C'...
        };
        
sio =  {'J25A','J33A','J35A','J36A',...
        'J43A','J44A','J57A','J65A','J73A',...
        'M01A','M02A','M04A','M05A','M07A','M08A',...
        'FS01B','FS14B',...
        'G02B','G10B','G12B','G27B','G28B','G37B',...
        'J09B','J10B','J18B','J20B',...
        'M11B','M12B','M14B',...
        'J25C','J33C','J52C','J53C',...
        'J61C','J65C','J68C','J73C',...
        'M01C','M02C','M03C','M04C','M05C','M07C','M08C',...
        };
        

trm = {'FN01A','FN05A','FN07A','FN08A',...
       'FN12A','FN14A','FN18A','J41A','J49A',...
       'FS03B','FS04B','FS08B',...
       'FS11B','FS12B','FS15B','FS17B','FS18B','FS19B',...
       'G09B','G25B','G33B',...
       'J17B','J25B','J33B',...
       'M09B','M10B','M13B','M18B',...
       'FN01C','FN02C','FN03C','FN04C','FN05C','FN06C','FN07C','FN08C','FN09C',...
       'FN10C','FN11C','FN12C','FN14C','FN17C','FN18C','FN19C',...
       'J41C','J49C','J57C',...
       };
   
switch sta
    case whoi
        IP = 'WHOI';
    case ldeo
        IP = 'LDEO';
    case sio
        IP = 'SIO';
end

if any(strcmp(sta,trm))
    note = 'TRM';
end

end