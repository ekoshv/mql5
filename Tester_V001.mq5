﻿//+------------------------------------------------------------------+
#property copyright "Copyright © 2023, Ehsan KhademOlama"
#property link      "https://www.lottosachen.de/"
#property version   "1.0"

#include <Trade\PositionInfo.mqh>
#include <Trade\Trade.mqh>
#include <Trade\SymbolInfo.mqh>
#include <Trade\AccountInfo.mqh>
#include <Trade\DealInfo.mqh>
#include <Trade\OrderInfo.mqh>
#include <Expert\Money\MoneyFixedMargin.mqh>
#include <Arrays\ArrayLong.mqh>
#include <Arrays\ArrayDouble.mqh>
#include <Arrays\ArrayInt.mqh>
#include <Math\Stat\Stat.mqh>
#include <mcarlo.mqh>
#include <CTradeStatistics.mqh>
//---
CPositionInfo  m_position;                   // object of CPositionInfo class
CTrade         m_trade;                      // object of CTrade class
CSymbolInfo    m_symbol;                     // object of CSymbolInfo class
CAccountInfo   m_account;                    // object of CAccountInfo class
CDealInfo      m_deal;                       // object of CDealInfo class
COrderInfo     m_order;                      // object of COrderInfo class
CMoneyFixedMargin *m_money;                  // object of CMoneyFixedMargin class
CTradeStatistics m_stat;                     // object of RealTime Trade Statistics
//+------------------------------------------------------------------+
//| Enum Lor or Risk                                                 |
//+------------------------------------------------------------------+
enum ENUM_LOT_OR_RISK
  {
   lots_min=0, // Lots Min
   lot=1,      // Constant lot
   risk=2,     // Fixed Margin
  };
enum ENUM_WEEKLY_OR_DAILY
  {
   daily = 0,
   weekly = 1,
  };
//+------------------------------------------------------------------+
//-----------------Money Withdrawal----------------------------------+
//+------------------------------------------------------------------+
input group "===Virtual Withdrawal Parameters===";//*******
static double         co_money    = 0;
static double         rem_money   = 0;
static ulong          wk          = 0;
sinput bool           MWE         = false;            // Money Withdrawal Enable
input ENUM_WEEKLY_OR_DAILY wdtype = daily;            // Type of Withdrawal (Daily or Weekly)
input double          wd_money    = 50;               // The Withdrawal Money: KxWM
input double          div         = 1.0;              // Extraction Fraction(0->1]
input double          tr_cof      = 1.0;              // Withdrawal Start/Initial 1->Inf
//+------------------------------------------------------------------+
//-----------------Global Parameters---------------------------------+
//+------------------------------------------------------------------+
input group "===Global Parameters===";//*******
input ulong   InpTakeProfit         = 500;        // Take-Profit, in points
input double  InpSLCoef             = 0.5;        // Stop-Loss Ratio (SL=TP*alp)
ulong   InpStopLoss           = (ulong)(InpTakeProfit*InpSLCoef);        // Stop Loss, in points
input double  InpTrCoef             = 0.5;        // Trailing Stop Ratio (Tr=SL*bet)
ulong   InpTrailingStop       = (ulong)(InpStopLoss*InpTrCoef);          // Trailing Stop, in points
input double  InpTrStCoef           = 0.5;        // Trailing Stop Step Ratio (TrSt=TP*gam)
ulong   InpTrailingStep       = (ulong)(InpTakeProfit*InpTrStCoef);      // Trailing Step, in points
input ulong   InpTrailingFrequency  = 10;         // Trailing, in seconds (< "10" -> only on a new bar)
input ENUM_LOT_OR_RISK InpLotOrRisk = lots_min;   // Money management: Lot OR Risk
input double   InpVolumeLotOrRisk   = 0.01;       // The value for "Money management"
input int      ttd                  = 2;          // Total trade in direction
input ulong    InpDeviation         = 10;         // Deviation, in points (1.00045-1.00055=10 points)
input bool     InpReverse           = false;      // Reverse Trading
input int      OKAcceptCount        = 4;          // Trade Trigger Active Counter
//input double   cutoffThresh         = 0.10;     // Counter CutOff Thresh (0.00->0.75)
input bool     InpCloseOpposite     = false;      // Close opposite
sinput bool    JustBuys             = true;       // Enable Buys
sinput bool    JustSells            = true;       // Enable Sells
sinput bool    enExit               = false;      // Enable Exit Strategies
input double   eThreshold           = 10;         // Exit Strategies Clousure threshold (Currency)
sinput bool   EnCloseDanger         = false;      // Close Dangerous Positions
input double   Dthresh              = 50;         // Close Loss more than
bool     InpPrintLog                = false;      // Print log
input ulong    InpMagic             = 17770101;     // magic number
ENUM_POSITION_TYPE TradeX  = InpReverse?POSITION_TYPE_SELL:POSITION_TYPE_BUY;
ENUM_POSITION_TYPE rTradeX = InpReverse?POSITION_TYPE_BUY:POSITION_TYPE_SELL;
//+------------------------------------------------------------------+
//--------------------Equity Extraction------------------------------|
//+------------------------------------------------------------------+
input group "===Equity Extraction Parameters===";//*******
sinput bool  EqEx                    = false; // Enable Equity Extraxtion
sinput bool  EqExPerFix              = true;  // Pertentage(True) or Fixed(Currency)
input double EqExValue               = 1;     // Extraction Value
double init_EqEx                     = 1;
//+------------------------------------------------------------------+
//--------------------Trade Expire-----------------------------------|
//+------------------------------------------------------------------+
input group "===Trade Expire Parameters===";//*******
sinput bool    ExpireEn             = true;       // Enable Position Expire
sinput bool    ExTermEn             = false;      // Enable Terminate Expired 2XTime
input  long    ExpireMinute         = 7/*Days*/*24*60;  // Expiration time (minutes)
input double   ExpireExtract        = 10;         // Profit Out (Currency)
input double   TermExtract          = 100;        // Termination Thresh (Currency)
//+------------------------------------------------------------------+
//--------------------Hedging----------------------------------------|
//+------------------------------------------------------------------+
input group "===Hedging Parameters===";//*******
sinput bool    HedgeEn                 = false;                   // Enable Hedging
sinput bool    ReHedgeEn               = false;                   // Reverse Hedge
sinput bool    Htype                   = true;                    // Static / Dynamic
input double   Inp_Hedge_DrawdownOpenX = 5;                       // Hedge Opening Drawdown (Points)
double         Hedge_DrawdownOpen      = Inp_Hedge_DrawdownOpenX;
input double   Inp_Hedge_Min_ProfitX   = 10;                      // Hedge Minimum Profit Extrction (Points)
double         Hedge_Min_Profit        = Inp_Hedge_Min_ProfitX;
input double   Hedge_Coef              = 1;                       // Hedge Volume Coef
CArrayLong     *HedgedPositionsID= new CArrayLong;
//+------------------------------------------------------------------+
//-----------------Shadow--------------------------------------------+
//+------------------------------------------------------------------+
input group "===Shadow Parameters===";//*******
sinput bool    shadowEn             = false;      // Shadow En
sinput bool    revShadow            = false;      // Reverse Shadow
input double   tpPercentX           = 0.5;        // Shadow Take Profit 0->1
input double   slPercentX           = 1.0;        // Shadow Stop Loss 0->1
double         tpPercent            = tpPercentX; // Shadow Take Profit 0->1
double         slPercent            = slPercentX; // Shadow Stop Loss 0->1
//+------------------------------------------------------------------+
//-----------------Fractal-------------------------------------------+
//+------------------------------------------------------------------+
input group "===Fractal Parameters===";//*******
input int    xcount              = 100;       //Fractal Counts
input double DelFrac             = 150;       //> Fractal Diff (Point)
input int FractFilt              =  10;       //Fractal Filter
input double fraka               = 0.01;         //Diff (%))
CArrayDouble     *FracsUpMA      = new CArrayDouble;
CArrayDouble     *FracsDownMA    = new CArrayDouble;

//+------------------------------------------------------------------+
//------------------------Brake Control------------------------------+
//+------------------------------------------------------------------+
input group "===Brake Control Parameters===";//*******
input bool              BrakeEn              =   false;    // Enable Brake
input bool              ActiveBrake          =   false;    // All/Lossy Positions
input double            BrakePerT            =   10.0;    // Brake (%From Max/Step)
input double            BrakeThreshold       =   20.0;    // Loss Brake Threshold (In Currency)
input double            brAlphaX             =    0.6;    // Brake Volume Coef (0.1->1)
double            SRref                =    2.5;    // Sharpe Ratio Ref
input ulong             gotoActiveCount      =      5;    // End Trading Brake Count
static double brAlpha;
//input double            BrakePerD            =   30.0;    // Dynamic Brake (From Max)
//input unsigned int      BrakeMinute          =   6*60;    // Minutes of Brake
static ulong brakeCountReal = 0;
static double BrVolume = 1;
static double SRXCoef = 1;
//+------------------------------------------------------------------+
//-----------------Time Control--------------------------------------+
//+------------------------------------------------------------------+
input group "===Time Control Parameters===";//*******
input int      st_hour           = 21;       // Trading Close Time (for Stocks 17-18)
input bool     InpTimeControl    = false;   // Use time control
input bool     InpTimeActive     = false;    // Active Method
input uchar    InpStartHour      = 8;       // Start hour
input uchar    InpEndHour        = 22;      // End hour
//+------------------------------------------------------------------+
//-----------------Super Optimizer Control--------------------------+
//+------------------------------------------------------------------+
input group "===Super Optimizer Parameters===";//*******
input double SharpeLow      = 1.2; // Sharpe Ratio Low Threshold
input double SharpeHighPlus = 1.5; // Sharpe Ratio High Threshold Plus
input double SRFac          = 3.0; // Sharpe Ratio Power
input double DWFactor       = 0.6; // Draw Down Factor
input double BlFac          =   4; // Balance Factor
input double BlLim          = 150; // Balance Thresh
input double EpRfFac        = 0.5; // Expected*Recovery Factor
input int    XtT            = 200; // Total Trade Minimum
input bool   starterEn      = false; // Enable Starter
input bool   NNtrainer      = false; // Enable NN
input bool   MaxBal         = false; // Enable Max Balance(pF*eP)
//+------------------------------------------------------------------+
//-----------------Monte Carlo Optimizer Control---------------------+
//+------------------------------------------------------------------+
input group "===Monte Carlo Optimizer Parameters===";//*******
input bool   MCOptEn           = false; // Monte Carlo Enable (Disable Super)
input int    noptpr            =     1; // optimization parameter variant, in [1,2,...,NOPTPRMAX(6)]
input double rmndmin           =   0.9; // drawdown restriction, in (0.0->1.0)
input double fwdsh             =   0.5; // share of deals in "future", in (0.0->1.0)
//+------------------------------------------------------------------+
//-----------------Internal------------------------------------------+
//+------------------------------------------------------------------+
double m_stop_loss               = 0.0;      // Stop Loss               -> double
double m_take_profit             = 0.0;      // Take Profit             -> double
double m_trailing_stop           = 0.0;      // Trailing Stop           -> double
double m_trailing_step           = 0.0;      // Trailing Step           -> double
int    MXEnBar = 3;
int    handle_iFractals;             // variable for storing the handle of the iFractals indicator

double   m_adjusted_point;                   // point value adjusted for 3 or 5 points
datetime m_last_trailing    = 0;             // "0" -> D'1970.01.01 00:00';
datetime m_prev_bars        = 0;             // "0" -> D'1970.01.01 00:00';
double ExtDrawdownOpen      = 0.0;           // Drawdown (opening a hedge) -> double
double ExtDrawdownClose     = 0.0;           // Drawdown (closing the hedge) -> double

const int    Nfr = 10;
double fractals_up[];
double fractals_down[];
int    fractals_up_number[];
int    fractals_down_number[];

static double p_MaxEq      = 1;
static double p_MaxBl      = 1;
static double p_Linear     = 1;
static double under_init   = 1;
static double under_MaxEq  = 1;
static double under_MaxBl  = 1;
static double under_Linear = 1;
static double init_Balance = 600;
static double init_Eq    = 600;
static double MaxEqNow = 1;
static double MaxBlNow = 1;
static double EqDiff = 1;
static double BlDiff = 1;
static double XDiff  = 1;
static ulong  dead_time = 0;
static ulong  sleep_time = 0;

static bool   InpOnlyOne  = false;
static ushort countBuys   = 0;
static ushort countSells  = 0;
static ulong  countClose  = 0;
static ulong  countProfitClose = 0;
static ulong  countXBuys  = 0;
static ulong  countXSells = 0;
static ulong  Positivised = 0;
static ulong  collectorCount = 0;
static ulong  brakeCount = 0;
static bool   OK_Buy = false;
static bool   OK_Sell = false;
static int    OK_Sell_Count = 0;
static int    OK_Buy_Count = 0;
static uint   CountTrades = 0;
uint          InitTradeCount = 5;
static int    BarCount = 0;
static double SRDiff = 1;
static double SRX = 0;
//--- the tactic is this: for positions we strictly monitor the result ***
//+------------------------------------------------------------------+
//| Structure Positions                                              |
//+------------------------------------------------------------------+
struct STRUCT_POSITION
  {
   ENUM_POSITION_TYPE pos_type;              // position type
   double            volume;                 // position volume (if "0.0" -> the lot is "Money management")
   double            lot_coefficient;        // lot coefficient
   bool              waiting_transaction;    // waiting transaction, "true" -> it's forbidden to trade, we expect a transaction
   ulong             waiting_order_ticket;   // waiting order ticket, ticket of the expected order
   bool              transaction_confirmed;  // transaction confirmed, "true" -> transaction confirmed
   string            comment;                // Order type, Normal or Hedge
   long              Identifier;
   double            takeProfitPercent;      // Level of Take Profit 0->1
   double            stopLossPercent;        // Level of Stop Loss 0->1
   //--- Constructor
                     STRUCT_POSITION()
     {
      pos_type                   = WRONG_VALUE;
      volume                     = 0.0;
      lot_coefficient            = 1.0;
      waiting_transaction        = false;
      waiting_order_ticket       = 0;
      transaction_confirmed      = false;
      comment                    = "Normal****";
      Identifier                 = -1;
      takeProfitPercent          = 1;
      stopLossPercent            = 1;
     }
  };
STRUCT_POSITION SPosition[];
//+------------------------------------------------------------------+
//| Expert initialization function                                   |
//+------------------------------------------------------------------+
int OnInit()
  {
   init_EqEx = m_account.Equity();
   m_stat.InitRiskFreeRate(0);
   brAlpha = (brAlphaX>=0.1)?brAlphaX:0.1;
   brAlpha = (brAlpha<=1)?brAlpha:1;
   init_Balance = m_account.Balance();
   init_Eq  = m_account.Equity();
   MaxEqNow = m_account.Equity();
   MaxBlNow = m_account.Balance();
   rem_money = init_Balance;
//---
   ArrayResize(fractals_up,Nfr);
   ArrayResize(fractals_down,Nfr);
   ArrayResize(fractals_up_number,Nfr);
   ArrayResize(fractals_down_number,Nfr);
   ArrayInitialize(fractals_up,EMPTY_VALUE);
   ArrayInitialize(fractals_down,EMPTY_VALUE);
   ArrayInitialize(fractals_up_number,-1);
   ArrayInitialize(fractals_down_number,-1);
//---
//---
   if(InpTrailingStop!=0 && InpTrailingStep==0)
     {
      string err_text="Trailing is not possible: parameter \"Trailing Step\" is zero!";
      //--- when testing, we will only output to the log about incorrect input parameters
      if(MQLInfoInteger(MQL_TESTER))
        {
         Print(", ERROR: ",err_text);
         //return(INIT_FAILED);
        }
      else // if the Expert Advisor is run on the chart, tell the user about the error
        {
         Alert(", ERROR: ",err_text);
         return(INIT_PARAMETERS_INCORRECT);
        }
     }
//---
   if(!m_symbol.Name(Symbol())) // sets symbol name
     {
      Print(", ERROR: CSymbolInfo.Name");
      //return(INIT_FAILED);
     }
   RefreshRates();
//---
   m_trade.SetExpertMagicNumber(InpMagic);
   m_trade.SetMarginMode();
   m_trade.SetTypeFillingBySymbol(m_symbol.Name());
   m_trade.SetDeviationInPoints(InpDeviation);
//--- tuning for 3 or 5 digits
   int digits_adjust=1;
   m_adjusted_point=m_symbol.Point()*digits_adjust;

   Print("!!!Digits: ",m_symbol.Digits());
   Print("!!!Point: ",m_symbol.Point());
   Print("!!!Adjust: ",m_adjusted_point);
   Print("!!!----------------------");

   m_stop_loss             = InpStopLoss              * m_adjusted_point;
   m_take_profit           = InpTakeProfit            * m_adjusted_point;
   m_trailing_stop         = InpTrailingStop          * m_adjusted_point;
   m_trailing_step         = InpTrailingStep          * m_adjusted_point;

   Hedge_DrawdownOpen      = Hedge_DrawdownOpen       * m_adjusted_point;
//Hedge_Max_Loss          = Hedge_Max_Loss           * m_adjusted_point;
   Hedge_Min_Profit        = Hedge_Min_Profit         * m_adjusted_point;
//--- check the input parameter "Lots"
   string err_text="";
   if(InpLotOrRisk==lot)
     {
      if(!CheckVolumeValue(InpVolumeLotOrRisk,err_text))
        {
         //--- when testing, we will only output to the log about incorrect input parameters
         if(MQLInfoInteger(MQL_TESTER))
           {
            Print(", ERROR: ",err_text);
            //return(INIT_FAILED);
           }
         else // if the Expert Advisor is run on the chart, tell the user about the error
           {
            Alert(", ERROR: ",err_text);
            return(INIT_PARAMETERS_INCORRECT);
           }
        }
     }
   else
      if(InpLotOrRisk==risk)
        {
         if(m_money!=NULL)
            delete m_money;
         m_money=new CMoneyFixedMargin;
         if(m_money!=NULL)
           {
            if(!m_money.Init(GetPointer(m_symbol),Period(),m_symbol.Point()*digits_adjust))
              {
               Print(", ERROR: CMoneySizeOptimized.Init");
               //return(INIT_FAILED);
              }
            m_money.Percent(InpVolumeLotOrRisk);
           }
         else
           {
            Print(", ERROR: Object CMoneySizeOptimized is NULL");
            //return(INIT_FAILED);
           }
        }
//---
//--- create handle of the indicator iFractals
   handle_iFractals=iFractals(m_symbol.Name(),Period());
//--- if the handle is not created
   if(handle_iFractals==INVALID_HANDLE)
     {
      //--- tell about the failure and output the error code
      PrintFormat("Failed to create handle of the iFractals indicator for the symbol %s/%s, error code %d",
                  m_symbol.Name(),
                  EnumToString(Period()),
                  GetLastError());
      //--- the indicator is stopped early
      //return(INIT_FAILED);
     }
//---
//---
   return(INIT_SUCCEEDED);
  }
//+------------------------------------------------------------------+
//| Expert deinitialization function                                 |
//+------------------------------------------------------------------+
void OnDeinit(const int reason)
  {
   Print("!!!Under Linear: ", under_Linear);
   Print("!!!Under Max Eq: ", under_MaxEq);
   Print("!!!Under Max Bl: ", under_MaxBl);
   Print("!!!Sleep Time: ",sleep_time);
   Print("!!!Dead Time: ",dead_time);
   Print("!!!Positivised: ",Positivised);
   Print("!!!Collected: ",collectorCount);
   Print("!!!Braked: ",brakeCountReal);
   Print("!!!Sharpe Difference: ", SRDiff);
   if(m_stat.Calculate())
     {
      PrintFormat("!!!Sharpe Ratio: %.2f",m_stat.SharpeRatio());
      PrintFormat("!!!LR Correlation: %.2f",m_stat.LRCorrelation());
      PrintFormat("!!!Z-Score: %.2f",m_stat.ZScore());
     }
   else
      Print(m_stat.GetLastErrorString());
//---
//---
   if(m_money!=NULL)
      delete m_money;
  }
//+------------------------------------------------------------------+
//| Expert tick function                                             |
//+------------------------------------------------------------------+
void OnTick()
  {

//---
   double Dlevel;
   if(FreezeStopsLevels(Dlevel) && EnCloseDanger)
     {
      ClosePositionsDangers(POSITION_TYPE_BUY,Dlevel,Dthresh);
      ClosePositionsDangers(POSITION_TYPE_SELL,Dlevel,Dthresh);
     }
//---
   datetime time_now = TimeCurrent();
   MqlDateTime date_now;
   TimeToStruct(time_now, date_now);

   if(date_now.year>2022 && date_now.mon>1)
      return;
//---
   if(EqEx)
     {
      double EqRat = 100*m_account.Profit()/init_EqEx;
      if(EqExPerFix)
        {
         if(EqRat>EqExValue)
           {
            double level;
            if(FreezeStopsLevels(level))
               ClosePositions(POSITION_TYPE_SELL,level);
            if(FreezeStopsLevels(level))
               ClosePositions(POSITION_TYPE_BUY,level);
            init_EqEx = m_account.Equity();
            Print("!!!Extraxted Per", init_EqEx);
           }//check value
        }//if perc of fix
      else
        {
         if(m_account.Profit()>EqExValue && m_account.Equity()>init_EqEx)
           {
            Print("!!!Extraxted Fix", m_account.Profit());
            double level;
            if(FreezeStopsLevels(level))
               ClosePositions(POSITION_TYPE_SELL,level);
            if(FreezeStopsLevels(level))
               ClosePositions(POSITION_TYPE_BUY,level);
            init_EqEx = m_account.Equity();
           }//check value
        }//else if perc or fix
     }//end of equity extraction
//---
   if(MWE)
     {
      static bool wd_en = true;
      switch(wdtype)
        {
         case daily:
            if(date_now.hour==st_hour && wd_en)
              {
               double alpm = m_account.Equity()-(tr_cof*init_Balance);
               if(alpm>0)
                 {
                  int kk = (int)floor((alpm*div/wd_money));
                  if(kk>0)
                    {

                     double level;
                     if(FreezeStopsLevels(level))
                        ClosePositions(POSITION_TYPE_SELL,level);
                     if(FreezeStopsLevels(level))
                        ClosePositions(POSITION_TYPE_BUY,level);

                     if(TesterWithdrawal(wd_money*kk))
                       {
                        wk++;
                        co_money += wd_money*kk;
                        rem_money = MathMax(m_account.Equity(),init_Balance);
                        under_MaxEq = 1;
                        //Print("!!!Rem: ",rem_money);
                        Print("!!!Money cashed Out today on ",
                              date_now.day_of_week," is ",wd_money*kk,
                              " & in total: ",co_money);
                        wd_en = false;
                       }
                    }
                 }
              }
            if(date_now.hour != st_hour)
               wd_en=true;
            break;
         //----------------------------------------------------------
         case weekly:
            if(date_now.day_of_week==5 && date_now.hour >= st_hour && wd_en)
              {
               double alpm = m_account.Equity()-(tr_cof*init_Balance);
               if(alpm>0)
                 {
                  int kk = (int)floor((alpm*div/wd_money));
                  if(kk>0)
                    {

                     double level;
                     if(FreezeStopsLevels(level))
                        ClosePositions(POSITION_TYPE_SELL,level);
                     if(FreezeStopsLevels(level))
                        ClosePositions(POSITION_TYPE_BUY,level);

                     if(TesterWithdrawal(wd_money*kk))
                       {
                        wk++;
                        co_money += wd_money*kk;
                        rem_money = MathMax(m_account.Equity(),init_Balance);
                        under_MaxEq = 1;
                        //Print("!!!Rem: ",rem_money);
                        Print("!!!Money cashed Out today on ",
                              date_now.day_of_week," is ",wd_money*kk,
                              " & in total: ",co_money);
                        wd_en = false;
                       }
                    }
                 }
              }
            if(date_now.day_of_week != 5)
               wd_en=true;
            break;
        }
     }
//--------------------------------------------------
//---
   if((date_now.day_of_week == 5 && date_now.hour >= st_hour)||
      (date_now.day_of_week == 6 || date_now.day_of_week == 0))
     {
      double level;
      if(FreezeStopsLevels(level))
         ClosePositionsProfitable(POSITION_TYPE_SELL,level,eThreshold);
      if(FreezeStopsLevels(level))
         ClosePositionsProfitable(POSITION_TYPE_BUY,level,eThreshold);
      return;
     }
//---
   if(InpTimeControl)
      if(!TimeControl())
        {
         if(InpTimeActive)
           {
            double level;
            if(FreezeStopsLevels(level))
               ClosePositions(POSITION_TYPE_SELL,level);
            if(FreezeStopsLevels(level))
               ClosePositions(POSITION_TYPE_BUY,level);
           }
         return;
        }
//---
   int      count_buys           = 0;
   double   volume_buys          = 0.0;
   double   volume_biggest_buys  = 0.0;
   int      count_sells          = 0;
   double   volume_sells         = 0.0;
   double   volume_biggest_sells = 0.0;
   int      count_hedge_buys     = 0;
   int      count_hedge_sells    = 0;
   int    count_X_buys     = 0;
   int    count_X_sells    = 0;

   CalculateAllPositions(count_buys,volume_buys,volume_biggest_buys,
                         count_sells,volume_sells,volume_biggest_sells,
                         count_hedge_buys,count_hedge_sells,
                         count_X_buys,
                         count_X_sells);
//---
   int size_need_position=ArraySize(SPosition);
   if(size_need_position>0)
     {
      for(int i=size_need_position-1; i>=0; i--)
        {
         if(SPosition[i].waiting_transaction)
           {
            if(!SPosition[i].transaction_confirmed)
              {
               if(InpPrintLog)
                  Print(", OK: ",
                        "transaction_confirmed: ",SPosition[i].transaction_confirmed);
               return;
              }
            else
               if(SPosition[i].transaction_confirmed)
                 {
                  ArrayRemove(SPosition,i,1);
                  return;
                 }
           }
         if(SPosition[i].pos_type==POSITION_TYPE_BUY)
           {
            if(InpCloseOpposite || InpOnlyOne)
              {
               if(InpCloseOpposite)
                 {
                  if(count_sells>0)
                    {
                     double level;
                     if(FreezeStopsLevels(level))
                        ClosePositions(POSITION_TYPE_SELL,level);
                     return;
                    }
                 }
               if(InpOnlyOne)
                 {
                  if(count_buys+count_sells==0)
                    {
                     double level;
                     if(FreezeStopsLevels(level))
                       {
                        SPosition[i].waiting_transaction=true;
                        OpenPosition(i,level);
                       }
                     return;
                    }
                  else
                     ArrayRemove(SPosition,i,1);
                  return;
                 }
              }
            double level;
            if(FreezeStopsLevels(level))
              {
               SPosition[i].waiting_transaction=true;
               OpenPosition(i,level);
              }
            return;
           }
         if(SPosition[i].pos_type==POSITION_TYPE_SELL)
           {
            if(InpCloseOpposite || InpOnlyOne)
              {

               CalculateAllPositions(count_buys,volume_buys,volume_biggest_buys,
                                     count_sells,volume_sells,volume_biggest_sells,
                                     count_hedge_buys,count_hedge_sells,
                                     count_X_buys,
                                     count_X_sells);
               if(InpCloseOpposite)
                 {
                  if(count_buys>0)
                    {
                     double level;
                     if(FreezeStopsLevels(level))
                        ClosePositions(POSITION_TYPE_BUY,level);
                     return;
                    }
                 }
               if(InpOnlyOne)
                 {
                  if(count_buys+count_sells==0)
                    {
                     double level;
                     if(FreezeStopsLevels(level))
                       {
                        SPosition[i].waiting_transaction=true;
                        OpenPosition(i,level);
                       }
                     return;
                    }
                  else
                     ArrayRemove(SPosition,i,1);
                  return;
                 }
              }
            double level;
            if(FreezeStopsLevels(level))
              {
               SPosition[i].waiting_transaction=true;
               OpenPosition(i,level);
              }
            return;
           }
        }
     }
//---

   if(InpTrailingFrequency>=10) // trailing no more than once every 10 seconds
     {
      datetime time_current=TimeCurrent();
      if(time_current-m_last_trailing>10)
        {
         double level;
         if(FreezeStopsLevels(level))
            Trailing(level);
         else
            return;
         m_last_trailing=time_current;
        }
     }


//-------------------Equity Check---------------------------------------------
   EqDiff = MaxEqNow>1?(m_account.Equity()-MaxEqNow)/MathSqrt(MaxEqNow):0;
   BlDiff = MaxBlNow>1?(m_account.Balance()-MaxBlNow)/MathSqrt(MaxBlNow):0;
   XDiff  = (MaxEqNow+MaxBlNow)>2?2*MathLog(MaxEqNow+MaxBlNow+1)*
            (m_account.Equity()-m_account.Balance())/
            MathSqrt(MaxEqNow+MaxBlNow):0;

   double EqP   = 100*m_account.Equity()/MaxEqNow;
   double BlP   = 100*m_account.Balance()/MaxBlNow;
   double EqInP = 100*m_account.Equity()/init_Eq;
   static bool BrakeOK = false;
   static bool BrakeFinished = false;
   static bool EndSession = false;
   static double PreMaxEq = init_Eq;
   static double PreMaxBl = init_Balance;

   if(BrakeEn)
     {
      if(MWE)
        {
         if(EqInP<(100-BrakePerT))
            BrakeOK = true;
        }
      else
        {
         if((EqP)<(100-BrakePerT) ||
            (BlP)<(100-BrakePerT))
           {
            //Print("!!!XEq: ",EqP," XBl: ",BlP);
            BrakeOK = true;
           }
        }

      if(BrakeFinished &&
         ((m_account.Equity()>PreMaxEq && m_account.Balance()>PreMaxBl)))
        {
         BrVolume = BrVolume/brAlpha;
         brakeCount--;
         BrakeFinished = false;
        }
     }

   if(BrakeOK)
     {
      if(ActiveBrake || EndSession)
        {
         EKOBraked();
         brakeCount++;
         brakeCountReal++;
         BrakeOK = false;
         BrakeFinished = true;
         return;
        }
      else
        {
         EKOBrakedX(BrakeThreshold);
         Print("!!!XEq: ",EqP," XBl: ",BlP);
         brakeCount++;
         brakeCountReal++;
         BrakeOK = false;
         PreMaxEq = MaxEqNow;
         MaxEqNow = m_account.Equity();
         PreMaxBl = MaxBlNow;
         MaxBlNow = m_account.Balance();
         BrVolume = brAlpha*BrVolume;
         BrakeFinished = true;
         if(brakeCount>=gotoActiveCount)
            EndSession = true;
         //return;
        }
      //return;
     }

   MaxEqNow = MathMax(MaxEqNow,m_account.Equity());
   MaxBlNow = MathMax(MaxBlNow,m_account.Balance());

//--- we work only at the time of the birth of new bar
   datetime time_0=iTime(m_symbol.Name(),Period(),0);
   if(time_0==m_prev_bars)
      return;
   m_prev_bars=time_0;

//--- Under Init Money Detection ---
   under_Linear += MathAbs(XDiff);

   if(EqDiff<0)
      under_MaxEq += MathAbs(EqDiff);

   if(BlDiff<0)
      under_MaxBl += MathAbs(BlDiff);

   double betm = (tr_cof*rem_money>0)?(m_account.Equity()-tr_cof*rem_money)/MathSqrt(tr_cof*rem_money):0;
   if(betm<0)
      under_init += MathAbs(betm);
//-----------------------------------------------------------------------------
//-------------------------------
   int      fcount_buys=0;
   double   fvolume_buys=0.0;
   double   fvolume_biggest_buys=0.0;
   int      fcount_sells=0;
   double   fvolume_sells=0.0;
   double   fvolume_biggest_sells=0.0;
   int      fcount_hedge_buys     = 0;
   int      fcount_hedge_sells    = 0;
   int      fcount_X_buys     = 0;
   int      fcount_X_sells    = 0;
   CalculateAllPositions(fcount_buys,fvolume_buys,fvolume_biggest_buys,
                         fcount_sells,fvolume_sells,fvolume_biggest_sells,
                         fcount_hedge_buys,fcount_hedge_sells,
                         fcount_X_buys, fcount_X_sells);

   if(fcount_buys+fcount_sells==0)
      dead_time++;
//-------------------------------
   double level;
   if(!FreezeStopsLevels(level))
     {
      m_prev_bars=0;
      return;
     }

   if(InpTrailingFrequency<10) // trailing only at the time of the birth of new bar
      Trailing(level);

//--- search for trading signals only at the time of the birth of new bar
   if(!SearchTradingSignals())
     {
      m_prev_bars=0;
      return;
     }
//---
   if(!EKOHedge())
     {
      m_prev_bars=0;
      return;
     }
//---
   if(ExpireEn)
     {
      double levelxp;
      if(FreezeStopsLevels(levelxp))
        {
         double X3 = m_account.Balance();
         EKOExpire(POSITION_TYPE_BUY,levelxp,ExpireExtract,ExpireMinute);
         if((m_account.Balance()-X3)>0)
           {
            countProfitClose++;
           }
        }
      if(FreezeStopsLevels(levelxp))
        {
         double X2 = m_account.Balance();
         EKOExpire(POSITION_TYPE_SELL,levelxp,ExpireExtract,ExpireMinute);
         if((m_account.Balance()-X2)>0)
           {
            countProfitClose++;
           }

        }
      if(ExTermEn)
        {
         double levelyp;
         if(FreezeStopsLevels(levelyp))
           {
            double X4 = m_account.Balance();
            EKOExpireTerm(POSITION_TYPE_BUY,levelyp,TermExtract,ExpireMinute);
            if((m_account.Balance()-X4)>0)
              {
               countProfitClose++;
              }

           }
         if(FreezeStopsLevels(levelyp))
           {
            double X5 = m_account.Balance();
            EKOExpireTerm(POSITION_TYPE_SELL,levelyp,TermExtract,ExpireMinute);
            if((m_account.Balance()-X5)>0)
              {
               countProfitClose++;
              }

           }
        }
     }
//---
  }
//+------------------------------------------------------------------+
//| TradeTransaction function                                        |
//+------------------------------------------------------------------+
void OnTradeTransaction(const MqlTradeTransaction &trans,
                        const MqlTradeRequest &request,
                        const MqlTradeResult &result)
  {
//--- get transaction type as enumeration value
   ENUM_TRADE_TRANSACTION_TYPE type=trans.type;
//--- if transaction is result of addition of the transaction in history
   if(type==TRADE_TRANSACTION_DEAL_ADD)
     {
      long     deal_ticket       =0;
      long     deal_order        =0;
      long     deal_time         =0;
      long     deal_time_msc     =0;
      long     deal_type         =-1;
      long     deal_entry        =-1;
      long     deal_magic        =0;
      long     deal_reason       =-1;
      long     deal_position_id  =0;
      double   deal_volume       =0.0;
      double   deal_price        =0.0;
      double   deal_commission   =0.0;
      double   deal_swap         =0.0;
      double   deal_profit       =0.0;
      string   deal_symbol       ="";
      string   deal_comment      ="";
      string   deal_external_id  ="";
      if(HistoryDealSelect(trans.deal))
        {
         deal_ticket       =HistoryDealGetInteger(trans.deal,DEAL_TICKET);
         deal_order        =HistoryDealGetInteger(trans.deal,DEAL_ORDER);
         deal_time         =HistoryDealGetInteger(trans.deal,DEAL_TIME);
         deal_time_msc     =HistoryDealGetInteger(trans.deal,DEAL_TIME_MSC);
         deal_type         =HistoryDealGetInteger(trans.deal,DEAL_TYPE);
         deal_entry        =HistoryDealGetInteger(trans.deal,DEAL_ENTRY);
         deal_magic        =HistoryDealGetInteger(trans.deal,DEAL_MAGIC);
         deal_reason       =HistoryDealGetInteger(trans.deal,DEAL_REASON);
         deal_position_id  =HistoryDealGetInteger(trans.deal,DEAL_POSITION_ID);

         deal_volume       =HistoryDealGetDouble(trans.deal,DEAL_VOLUME);
         deal_price        =HistoryDealGetDouble(trans.deal,DEAL_PRICE);
         deal_commission   =HistoryDealGetDouble(trans.deal,DEAL_COMMISSION);
         deal_swap         =HistoryDealGetDouble(trans.deal,DEAL_SWAP);
         deal_profit       =HistoryDealGetDouble(trans.deal,DEAL_PROFIT);

         deal_symbol       =HistoryDealGetString(trans.deal,DEAL_SYMBOL);
         deal_comment      =HistoryDealGetString(trans.deal,DEAL_COMMENT);
         deal_external_id  =HistoryDealGetString(trans.deal,DEAL_EXTERNAL_ID);
        }
      else
         return;
      ENUM_DEAL_ENTRY enum_deal_entry=(ENUM_DEAL_ENTRY)deal_entry;
      if(deal_symbol==m_symbol.Name() && deal_magic==InpMagic)
        {
         if(deal_type==DEAL_TYPE_BUY || deal_type==DEAL_TYPE_SELL)
           {
            int size_need_position=ArraySize(SPosition);
            if(size_need_position>0)
              {
               for(int i=0; i<size_need_position; i++)
                 {
                  if(SPosition[i].waiting_transaction)
                     if(SPosition[i].waiting_order_ticket==deal_order)
                       {
                        Print(__FUNCTION__," Transaction confirmed");
                        SPosition[i].transaction_confirmed=true;
                        break;
                       }
                 }
              }
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| Refreshes the symbol quotes data                                 |
//+------------------------------------------------------------------+
bool RefreshRates()
  {
//--- refresh rates
   if(!m_symbol.RefreshRates())
     {
      if(InpPrintLog)
         Print(", ERROR: ","RefreshRates error");
      return(false);
     }
//--- protection against the return value of "zero"
   if(m_symbol.Ask()==0 || m_symbol.Bid()==0)
     {
      if(InpPrintLog)
         Print(", ERROR: ","Ask == 0.0 OR Bid == 0.0");
      return(false);
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| Check the correctness of the position volume                     |
//+------------------------------------------------------------------+
bool CheckVolumeValue(double volume,string &error_description)
  {
//--- minimal allowed volume for trade operations
   double min_volume=m_symbol.LotsMin();
   if(volume<min_volume)
     {
      error_description=StringFormat("Volume is less than the minimal allowed SYMBOL_VOLUME_MIN=%.2f",min_volume);
      return(false);
     }
//--- maximal allowed volume of trade operations
   double max_volume=m_symbol.LotsMax();
   if(volume>max_volume)
     {
      error_description=StringFormat("Volume is greater than the maximal allowed SYMBOL_VOLUME_MAX=%.2f",max_volume);
      return(false);
     }
//--- get minimal step of volume changing
   double volume_step=m_symbol.LotsStep();
   int ratio=(int)MathRound(volume/volume_step);
   if(MathAbs(ratio*volume_step-volume)>0.0000001)
     {
      error_description=StringFormat("Volume is not a multiple of the minimal step SYMBOL_VOLUME_STEP=%.2f, the closest correct volume is %.2f",
                                     volume_step,ratio*volume_step);
      return(false);
     }
   error_description="Correct volume value";
   return(true);
  }
//+------------------------------------------------------------------+
//| Lot Check                                                        |
//+------------------------------------------------------------------+
double LotCheck(double lots,CSymbolInfo &symbol)
  {
//--- calculate maximum volume
   double volume=NormalizeDouble(lots,2);
   double stepvol=symbol.LotsStep();
   if(stepvol>0.0)
      volume=stepvol*MathFloor(volume/stepvol);
//---
   double minvol=symbol.LotsMin();
   if(volume<minvol)
      volume=0.0;
//---
   double maxvol=symbol.LotsMax();
   if(volume>maxvol)
      volume=maxvol;
   return(volume);
  }
//+------------------------------------------------------------------+
//| Check Freeze and Stops levels                                    |
//+------------------------------------------------------------------+
bool FreezeStopsLevels(double &level)
  {
//--- check Freeze and Stops levels
   /*
      Type of order/position  |  Activation price  |  Check
      ------------------------|--------------------|--------------------------------------------
      Buy Limit order         |  Ask               |  Ask-OpenPrice  >= SYMBOL_TRADE_FREEZE_LEVEL
      Buy Stop order          |  Ask               |  OpenPrice-Ask  >= SYMBOL_TRADE_FREEZE_LEVEL
      Sell Limit order        |  Bid               |  OpenPrice-Bid  >= SYMBOL_TRADE_FREEZE_LEVEL
      Sell Stop order         |  Bid               |  Bid-OpenPrice  >= SYMBOL_TRADE_FREEZE_LEVEL
      Buy position            |  Bid               |  TakeProfit-Bid >= SYMBOL_TRADE_FREEZE_LEVEL
                              |                    |  Bid-StopLoss   >= SYMBOL_TRADE_FREEZE_LEVEL
      Sell position           |  Ask               |  Ask-TakeProfit >= SYMBOL_TRADE_FREEZE_LEVEL
                              |                    |  StopLoss-Ask   >= SYMBOL_TRADE_FREEZE_LEVEL

      Buying is done at the Ask price                 |  Selling is done at the Bid price
      ------------------------------------------------|----------------------------------
      TakeProfit        >= Bid                        |  TakeProfit        <= Ask
      StopLoss          <= Bid                        |  StopLoss          >= Ask
      TakeProfit - Bid  >= SYMBOL_TRADE_STOPS_LEVEL   |  Ask - TakeProfit  >= SYMBOL_TRADE_STOPS_LEVEL
      Bid - StopLoss    >= SYMBOL_TRADE_STOPS_LEVEL   |  StopLoss - Ask    >= SYMBOL_TRADE_STOPS_LEVEL
   */
   if(!RefreshRates() || !m_symbol.Refresh())
      return(false);
//--- FreezeLevel -> for pending order and modification
   double freeze_level=m_symbol.FreezeLevel()*m_symbol.Point();
   if(freeze_level==0.0)
      freeze_level=(m_symbol.Ask()-m_symbol.Bid())*3.0;
//--- StopsLevel -> for TakeProfit and StopLoss
   double stop_level=m_symbol.StopsLevel()*m_symbol.Point();
   if(stop_level==0.0)
      stop_level=(m_symbol.Ask()-m_symbol.Bid())*3.0;

   if(freeze_level<=0.0 || stop_level<=0.0)
      return(false);

   level=(freeze_level>stop_level)?freeze_level:stop_level;

   double spread=m_symbol.Spread()*m_symbol.Point()*3.0;
   level=(level>spread)?level:spread;
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| Open position                                                    |
//|   double stop_loss                                               |
//|      -> pips * m_adjusted_point (if "0.0" -> the m_stop_loss)    |
//|   double take_profit                                             |
//|      -> pips * m_adjusted_point (if "0.0" -> the m_take_profit)  |
//+------------------------------------------------------------------+
void OpenPosition(const int index,const double level)
  {
   /*
      Buying is done at the Ask price                 |  Selling is done at the Bid price
      ------------------------------------------------|----------------------------------
      TakeProfit        >= Bid                        |  TakeProfit        <= Ask
      StopLoss          <= Bid                        |  StopLoss          >= Ask
      TakeProfit - Bid  >= SYMBOL_TRADE_STOPS_LEVEL   |  Ask - TakeProfit  >= SYMBOL_TRADE_STOPS_LEVEL
      Bid - StopLoss    >= SYMBOL_TRADE_STOPS_LEVEL   |  StopLoss - Ask    >= SYMBOL_TRADE_STOPS_LEVEL
   */
//--- buy
   if(SPosition[index].pos_type==POSITION_TYPE_BUY)
     {
      double sl = -10;
      double alp = 1;
      while(sl<0)
        {
         sl=(m_stop_loss==0.0)?0.0:m_symbol.Ask()-(alp*m_stop_loss*SPosition[index].stopLossPercent);
         if(sl>0.0)
            if(m_symbol.Bid()-sl<level)
               sl=m_symbol.Bid()-level;
         alp = 3*alp/4;
         //Print("!!!Positivized SL!!!");
         Positivised++;
        }
      double tp=(m_take_profit==0.0)?0.0:m_symbol.Ask()+(m_take_profit*SPosition[index].takeProfitPercent);
      if(tp>0.0)
         if(tp-m_symbol.Ask()<level)
            tp=m_symbol.Ask()+level;

      OpenBuy(index,sl,tp);
      return;
     }
//--- sell
   if(SPosition[index].pos_type==POSITION_TYPE_SELL)
     {
      double sl=(m_stop_loss==0.0)?0.0:m_symbol.Bid()+(m_stop_loss*SPosition[index].stopLossPercent);
      if(sl>0.0)
         if(sl-m_symbol.Ask()<level)
            sl=m_symbol.Ask()+level;

      double tp = -10;
      double bet = 1;
      while(tp<0)
        {
         tp=(m_take_profit==0.0)?0.0:m_symbol.Bid()-(bet*m_take_profit*SPosition[index].takeProfitPercent);
         if(tp>0.0)
            if(m_symbol.Bid()-tp<level)
               tp=m_symbol.Bid()-level;
         bet = 3*bet/4;
         //Print("!!!Positivized TP!!!");
         Positivised++;
        }
      OpenSell(index,sl,tp);
      return;
     }
  }
//+------------------------------------------------------------------+
//| Open Buy position                                                |
//+------------------------------------------------------------------+
void OpenBuy(const int index,double sl,double tp)
  {
   sl=m_symbol.NormalizePrice(sl);
   tp=m_symbol.NormalizePrice(tp);

   double long_lot=0.0;
   if(SPosition[index].volume>0.0)
      long_lot=SPosition[index].volume;
   else
     {
      if(InpLotOrRisk==risk)
        {
         long_lot=m_money.CheckOpenLong(m_symbol.Ask(),sl);
         if(InpPrintLog)
            Print(", OK: ","sl=",DoubleToString(sl,m_symbol.Digits()),
                  ", CheckOpenLong: ",DoubleToString(long_lot,2),
                  ", Balance: ",    DoubleToString(m_account.Balance(),2),
                  ", Equity: ",     DoubleToString(m_account.Equity(),2),
                  ", FreeMargin: ", DoubleToString(m_account.FreeMargin(),2));
         if(long_lot==0.0)
           {
            ArrayRemove(SPosition,index,1);
            if(InpPrintLog)
               Print(", ERROR: ","CMoneyFixedMargin.CheckOpenLong returned the value of 0.0");
            return;
           }
        }
      else
         if(InpLotOrRisk==lot)
            long_lot=InpVolumeLotOrRisk;
         else
            if(InpLotOrRisk==lots_min)
               long_lot=m_symbol.LotsMin();
            else
              {
               ArrayRemove(SPosition,index,1);
               return;
              }
     }
   if(SPosition[index].lot_coefficient>0.0)
     {
      long_lot=LotCheck(long_lot*SPosition[index].lot_coefficient,
                        m_symbol);
      if(long_lot==0)
        {
         ArrayRemove(SPosition,index,1);
         if(InpPrintLog)
            Print(", ERROR: ","LotCheck returned the 0.0");
         return;
        }
     }
   if(m_symbol.LotsLimit()>0.0)
     {
      int      count_buys           = 0;
      double   volume_buys          = 0.0;
      double   volume_biggest_buys  = 0.0;
      int      count_sells          = 0;
      double   volume_sells         = 0.0;
      double   volume_biggest_sells = 0.0;
      int      count_hedge_buys    = 0;
      int      count_hedge_sells   = 0;
      int    count_X_buys     = 0;
      int    count_X_sells    = 0;
      CalculateAllPositions(count_buys,volume_buys,volume_biggest_buys,
                            count_sells,volume_sells,volume_biggest_sells,
                            count_hedge_buys,count_hedge_sells,
                            count_X_buys,count_X_sells);

      if(volume_buys+volume_sells+long_lot>m_symbol.LotsLimit())
        {
         ArrayRemove(SPosition,index,1);
         if(InpPrintLog)
            Print(", ERROR: ","#0 Buy, Volume Buy (",DoubleToString(volume_buys,2),
                  ") + Volume Sell (",DoubleToString(volume_sells,2),
                  ") + Volume long (",DoubleToString(long_lot,2),
                  ") > Lots Limit (",DoubleToString(m_symbol.LotsLimit(),2),")");
         return;
        }
     }
//--- check volume before OrderSend to avoid "not enough money" error (CTrade)
   double free_margin_check=m_account.FreeMarginCheck(m_symbol.Name(),
                            ORDER_TYPE_BUY,
                            long_lot,
                            m_symbol.Ask());
   double margin_check=m_account.MarginCheck(m_symbol.Name(),
                       ORDER_TYPE_BUY,
                       long_lot,
                       m_symbol.Ask());

   if(free_margin_check>margin_check)
     {
      if(m_trade.Buy(long_lot,m_symbol.Name(),
                     m_symbol.Ask(),
                     sl,tp,
                     SPosition[index].comment)) // CTrade::Buy -> "true"
        {
         //if(SPosition[index].comment == "Hedge_Buy" || SPosition[index].comment == "Hedge_Sell")
         //   HedgedPositionsID.Add(SPosition[index].Identifier);
         if(m_trade.ResultDeal()==0)
           {
            if(m_trade.ResultRetcode()==10009) // trade order went to the exchange
              {
               SPosition[index].waiting_transaction=true;
               SPosition[index].waiting_order_ticket=m_trade.ResultOrder();
              }
            else
              {
               SPosition[index].waiting_transaction=false;
               if(InpPrintLog)
                  Print(", ERROR: ","#1 Buy -> false. Result Retcode: ",m_trade.ResultRetcode(),
                        ", description of result: ",m_trade.ResultRetcodeDescription());
              }
            if(InpPrintLog)
               PrintResultTrade(m_trade,m_symbol);
           }
         else
           {
            if(m_trade.ResultRetcode()==10009)
              {
               SPosition[index].waiting_transaction=true;
               SPosition[index].waiting_order_ticket=m_trade.ResultOrder();
              }
            else
              {
               SPosition[index].waiting_transaction=false;
               if(InpPrintLog)
                  Print(", OK: ","#2 Buy -> true. Result Retcode: ",m_trade.ResultRetcode(),
                        ", description of result: ",m_trade.ResultRetcodeDescription());
              }
            if(InpPrintLog)
               PrintResultTrade(m_trade,m_symbol);
           }
        }
      else
        {
         SPosition[index].waiting_transaction=false;
         if(InpPrintLog)
            Print(", ERROR: ","#3 Buy -> false. Result Retcode: ",m_trade.ResultRetcode(),
                  ", description of result: ",m_trade.ResultRetcodeDescription());
         if(InpPrintLog)
            PrintResultTrade(m_trade,m_symbol);
        }
     }
   else
     {
      ArrayRemove(SPosition,index,1);
      if(InpPrintLog)
         Print(", ERROR: ","CAccountInfo.FreeMarginCheck returned the value ",DoubleToString(free_margin_check,2));
      return;
     }
//---
  }
//+------------------------------------------------------------------+
//| Open Sell position                                               |
//+------------------------------------------------------------------+
void OpenSell(const int index,double sl,double tp)
  {
   sl=m_symbol.NormalizePrice(sl);
   tp=m_symbol.NormalizePrice(tp);

   double short_lot=0.0;
   if(SPosition[index].volume>0.0)
      short_lot=SPosition[index].volume;
   else
     {
      if(InpLotOrRisk==risk)
        {
         short_lot=m_money.CheckOpenShort(m_symbol.Bid(),sl);
         if(InpPrintLog)
            Print(", OK: ","sl=",DoubleToString(sl,m_symbol.Digits()),
                  ", CheckOpenLong: ",DoubleToString(short_lot,2),
                  ", Balance: ",    DoubleToString(m_account.Balance(),2),
                  ", Equity: ",     DoubleToString(m_account.Equity(),2),
                  ", FreeMargin: ", DoubleToString(m_account.FreeMargin(),2));
         if(short_lot==0.0)
           {
            ArrayRemove(SPosition,index,1);
            if(InpPrintLog)
               Print(", ERROR: ","CMoneyFixedMargin.CheckOpenShort returned the value of \"0.0\"");
            return;
           }
        }
      else
         if(InpLotOrRisk==lot)
            short_lot=InpVolumeLotOrRisk;
         else
            if(InpLotOrRisk==lots_min)
               short_lot=m_symbol.LotsMin();
            else
              {
               ArrayRemove(SPosition,index,1);
               return;
              }
     }
   if(SPosition[index].lot_coefficient>0.0)
     {
      short_lot=LotCheck(short_lot*SPosition[index].lot_coefficient,m_symbol);
      if(short_lot==0)
        {
         ArrayRemove(SPosition,index,1);
         if(InpPrintLog)
            Print(", ERROR: ","LotCheck returned the 0.0");
         return;
        }
     }
   if(m_symbol.LotsLimit()>0.0)
     {
      int      count_buys           = 0;
      double   volume_buys          = 0.0;
      double   volume_biggest_buys  = 0.0;
      int      count_sells          = 0;
      double   volume_sells         = 0.0;
      double   volume_biggest_sells = 0.0;
      int      count_hedge_buys     = 0;
      int      count_hedge_sells    = 0;
      int    count_X_buys     = 0;
      int    count_X_sells    = 0;
      CalculateAllPositions(count_buys,volume_buys,volume_biggest_buys,
                            count_sells,volume_sells,volume_biggest_sells,
                            count_hedge_buys,count_hedge_sells,
                            count_X_buys, count_X_sells);
      if(volume_buys+volume_sells+short_lot>m_symbol.LotsLimit())
        {
         ArrayRemove(SPosition,index,1);
         if(InpPrintLog)
            Print(", ERROR: ","#0 Buy, Volume Buy (",DoubleToString(volume_buys,2),
                  ") + Volume Sell (",DoubleToString(volume_sells,2),
                  ") + Volume short (",DoubleToString(short_lot,2),
                  ") > Lots Limit (",DoubleToString(m_symbol.LotsLimit(),2),")");
         return;
        }
     }
//--- check volume before OrderSend to avoid "not enough money" error (CTrade)
   double free_margin_check=m_account.FreeMarginCheck(m_symbol.Name(),
                            ORDER_TYPE_SELL,
                            short_lot,
                            m_symbol.Bid());
   double margin_check=m_account.MarginCheck(m_symbol.Name(),
                       ORDER_TYPE_SELL,
                       short_lot,
                       m_symbol.Bid());
   if(free_margin_check>margin_check)
     {
      if(m_trade.Sell(short_lot,m_symbol.Name(),
                      m_symbol.Bid(),
                      sl,tp,
                      SPosition[index].comment)) // CTrade::Sell -> "true"
        {
         //if(SPosition[index].comment == "Hedge_Buy" || SPosition[index].comment == "Hedge_Sell")
         //   HedgedPositionsID.Add(SPosition[index].Identifier);
         //Print("!!!sl: ",sl);
         if(m_trade.ResultDeal()==0)
           {
            if(m_trade.ResultRetcode()==10009) // trade order went to the exchange
              {
               SPosition[index].waiting_transaction=true;
               SPosition[index].waiting_order_ticket=m_trade.ResultOrder();
              }
            else
              {
               SPosition[index].waiting_transaction=false;
               if(InpPrintLog)
                  Print(", ERROR: ","#1 Sell -> false. Result Retcode: ",m_trade.ResultRetcode(),
                        ", description of result: ",m_trade.ResultRetcodeDescription());
              }
            if(InpPrintLog)
               PrintResultTrade(m_trade,m_symbol);
           }
         else
           {
            if(m_trade.ResultRetcode()==10009)
              {
               SPosition[index].waiting_transaction=true;
               SPosition[index].waiting_order_ticket=m_trade.ResultOrder();
              }
            else
              {
               SPosition[index].waiting_transaction=false;
               if(InpPrintLog)
                  Print(", OK: ","#2 Sell -> true. Result Retcode: ",m_trade.ResultRetcode(),
                        ", description of result: ",m_trade.ResultRetcodeDescription());
              }
            if(InpPrintLog)
               PrintResultTrade(m_trade,m_symbol);
           }
        }
      else
        {
         SPosition[index].waiting_transaction=false;
         if(InpPrintLog)
            Print(", ERROR: ","#3 Sell -> false. Result Retcode: ",m_trade.ResultRetcode(),
                  ", description of result: ",m_trade.ResultRetcodeDescription());
         if(InpPrintLog)
            PrintResultTrade(m_trade,m_symbol);
        }
     }
   else
     {
      ArrayRemove(SPosition,index,1);
      if(InpPrintLog)
         Print(", ERROR: ","CAccountInfo.FreeMarginCheck returned the value ",DoubleToString(free_margin_check,2));
      return;
     }
//---
  }
//+------------------------------------------------------------------+
//| Print CTrade result                                              |
//+------------------------------------------------------------------+
void PrintResultTrade(CTrade &trade,CSymbolInfo &symbol)
  {
   Print(", Symbol: ",symbol.Name()+", "+
         "Code of request result: "+IntegerToString(trade.ResultRetcode())+", "+
         "Code of request result as a string: "+trade.ResultRetcodeDescription());
   Print("Deal ticket: "+IntegerToString(trade.ResultDeal())+", "+
         "Order ticket: "+IntegerToString(trade.ResultOrder())+", "+
         "Order retcode external: "+IntegerToString(trade.ResultRetcodeExternal())+", "+
         "Volume of deal or order: "+DoubleToString(trade.ResultVolume(),2));
   Print("Price, confirmed by broker: "+DoubleToString(trade.ResultPrice(),symbol.Digits())+", "+
         "Current bid price: "+DoubleToString(symbol.Bid(),symbol.Digits())+" (the requote): "+DoubleToString(trade.ResultBid(),symbol.Digits())+", "+
         "Current ask price: "+DoubleToString(symbol.Ask(),symbol.Digits())+" (the requote): "+DoubleToString(trade.ResultAsk(),symbol.Digits()));
   Print("Broker comment: "+trade.ResultComment());
  }
//+------------------------------------------------------------------+
//| Get value of buffers                                             |
//+------------------------------------------------------------------+
bool iGetArray(const int handle,const int buffer,const int start_pos,
               const int count,double &arr_buffer[])
  {
   bool result=true;
   if(!ArrayIsDynamic(arr_buffer))
     {
      if(InpPrintLog)
         PrintFormat("ERROR! EA: %s, FUNCTION: %s, this a no dynamic array!",__FILE__,__FUNCTION__);
      return(false);
     }
   ArrayFree(arr_buffer);
//--- reset error code
   ResetLastError();
//--- fill a part of the iBands array with values from the indicator buffer
   int copied=CopyBuffer(handle,buffer,start_pos,count,arr_buffer);
   if(copied!=count)
     {
      //--- if the copying fails, tell the error code
      if(InpPrintLog)
         PrintFormat("ERROR! EA: %s, FUNCTION: %s, amount to copy: %d, copied: %d, error code %d",
                     __FILE__,__FUNCTION__,count,copied,GetLastError());
      //--- quit with zero result - it means that the indicator is considered as not calculated
      return(false);
     }
   return(result);
  }
//+------------------------------------------------------------------+
//| Trailing                                                         |
//|   InpTrailingStop: min distance from price to Stop Loss          |
//+------------------------------------------------------------------+
void Trailing(const double stop_level)
  {
   /*
      Buying is done at the Ask price                 |  Selling is done at the Bid price
      ------------------------------------------------|----------------------------------
      TakeProfit        >= Bid                        |  TakeProfit        <= Ask
      StopLoss          <= Bid                        |  StopLoss          >= Ask
      TakeProfit - Bid  >= SYMBOL_TRADE_STOPS_LEVEL   |  Ask - TakeProfit  >= SYMBOL_TRADE_STOPS_LEVEL
      Bid - StopLoss    >= SYMBOL_TRADE_STOPS_LEVEL   |  StopLoss - Ask    >= SYMBOL_TRADE_STOPS_LEVEL
   */
   if(InpTrailingStop==0)
      return;
   for(int i=PositionsTotal()-1; i>=0; i--) // returns the number of open positions
      if(m_position.SelectByIndex(i))
         if(m_position.Symbol()==m_symbol.Name() && m_position.Magic()==InpMagic && m_position.Comment()=="Normal****")
           {
            double price_current = m_position.PriceCurrent();
            double price_open    = m_position.PriceOpen();
            double stop_loss     = m_position.StopLoss();
            double take_profit   = m_position.TakeProfit();
            double ask           = m_symbol.Ask();
            double bid           = m_symbol.Bid();
            //---
            if(m_position.PositionType()==POSITION_TYPE_BUY)
              {
               if(price_current-price_open>m_trailing_stop+m_trailing_step)
                  if(stop_loss<price_current-(m_trailing_stop+m_trailing_step))
                     if(m_trailing_stop>=stop_level && (take_profit-bid>=stop_level || take_profit==0.0))
                       {
                        if(!m_trade.PositionModify(m_position.Ticket(),
                                                   m_symbol.NormalizePrice(price_current-m_trailing_stop),
                                                   take_profit))
                           if(InpPrintLog)
                              Print(", ERROR: ","Modify BUY ",m_position.Ticket(),
                                    " Position -> false. Result Retcode: ",m_trade.ResultRetcode(),
                                    ", description of result: ",m_trade.ResultRetcodeDescription());
                        if(InpPrintLog)
                          {
                           RefreshRates();
                           m_position.SelectByIndex(i);
                           PrintResultModify(m_trade,m_symbol,m_position);
                          }
                        continue;
                       }
              }
            else
              {
               if(price_open-price_current>m_trailing_stop+m_trailing_step)
                  if((stop_loss>(price_current+(m_trailing_stop+m_trailing_step))) || (stop_loss==0))
                     if(m_trailing_stop>=stop_level && ask-take_profit>=stop_level)
                       {
                        if(!m_trade.PositionModify(m_position.Ticket(),
                                                   m_symbol.NormalizePrice(price_current+m_trailing_stop),
                                                   take_profit))
                           if(InpPrintLog)
                              Print(", ERROR: ","Modify SELL ",m_position.Ticket(),
                                    " Position -> false. Result Retcode: ",m_trade.ResultRetcode(),
                                    ", description of result: ",m_trade.ResultRetcodeDescription());
                        if(InpPrintLog)
                          {
                           RefreshRates();
                           m_position.SelectByIndex(i);
                           PrintResultModify(m_trade,m_symbol,m_position);
                          }
                       }
              }
           }
  }
//+------------------------------------------------------------------+
//| Print CTrade result                                              |
//+------------------------------------------------------------------+
void PrintResultModify(CTrade &trade,CSymbolInfo &symbol,CPositionInfo &position)
  {
//Print("File: ",__FILE__,", symbol: ",symbol.Name());
   Print("Code of request result: "+IntegerToString(trade.ResultRetcode()));
   Print("code of request result as a string: "+trade.ResultRetcodeDescription());
   Print("Deal ticket: "+IntegerToString(trade.ResultDeal()));
   Print("Order ticket: "+IntegerToString(trade.ResultOrder()));
   Print("Volume of deal or order: "+DoubleToString(trade.ResultVolume(),2));
   Print("Price, confirmed by broker: "+DoubleToString(trade.ResultPrice(),symbol.Digits()));
   Print("Current bid price: "+DoubleToString(symbol.Bid(),symbol.Digits())+" (the requote): "+DoubleToString(trade.ResultBid(),symbol.Digits()));
   Print("Current ask price: "+DoubleToString(symbol.Ask(),symbol.Digits())+" (the requote): "+DoubleToString(trade.ResultAsk(),symbol.Digits()));
   Print("Broker comment: "+trade.ResultComment());
   Print("Freeze Level: "+DoubleToString(symbol.FreezeLevel(),0),", Stops Level: "+DoubleToString(symbol.StopsLevel(),0));
   Print("Price of position opening: "+DoubleToString(position.PriceOpen(),symbol.Digits()));
   Print("Price of position's Stop Loss: "+DoubleToString(position.StopLoss(),symbol.Digits()));
   Print("Price of position's Take Profit: "+DoubleToString(position.TakeProfit(),symbol.Digits()));
   Print("Current price by position: "+DoubleToString(position.PriceCurrent(),symbol.Digits()));
  }
//+------------------------------------------------------------------+
//| Close positions                                                  |
//+------------------------------------------------------------------+
void ClosePositions(const ENUM_POSITION_TYPE pos_type,const double level)
  {
   /*
      Buying is done at the Ask price                 |  Selling is done at the Bid price
      ------------------------------------------------|----------------------------------
      TakeProfit        >= Bid                        |  TakeProfit        <= Ask
      StopLoss          <= Bid                        |  StopLoss          >= Ask
      TakeProfit - Bid  >= SYMBOL_TRADE_STOPS_LEVEL   |  Ask - TakeProfit  >= SYMBOL_TRADE_STOPS_LEVEL
      Bid - StopLoss    >= SYMBOL_TRADE_STOPS_LEVEL   |  StopLoss - Ask    >= SYMBOL_TRADE_STOPS_LEVEL
   */
   for(int i=PositionsTotal()-1; i>=0; i--) // returns the number of current positions
      if(m_position.SelectByIndex(i)) // selects the position by index for further access to its properties
         if(m_position.Symbol()==m_symbol.Name() && m_position.Magic()==InpMagic)
            if(m_position.PositionType()==pos_type)
              {
               if(m_position.PositionType()==POSITION_TYPE_BUY)
                 {
                  bool take_profit_level=(m_position.TakeProfit()!=0.0 && m_position.TakeProfit()-m_position.PriceCurrent()>=level) || m_position.TakeProfit()==0.0;
                  bool stop_loss_level=(m_position.StopLoss()!=0.0 && m_position.PriceCurrent()-m_position.StopLoss()>=level) || m_position.StopLoss()==0.0;
                  if(take_profit_level && stop_loss_level)
                     if(!m_trade.PositionClose(m_position.Ticket())) // close a position by the specified m_symbol
                        Print(", ERROR: ","CTrade.PositionClose ",m_position.Ticket());
                 }
               if(m_position.PositionType()==POSITION_TYPE_SELL)
                 {
                  bool take_profit_level=(m_position.TakeProfit()!=0.0 && m_position.PriceCurrent()-m_position.TakeProfit()>=level) || m_position.TakeProfit()==0.0;
                  bool stop_loss_level=(m_position.StopLoss()!=0.0 && m_position.StopLoss()-m_position.PriceCurrent()>=level) || m_position.StopLoss()==0.0;
                  if(take_profit_level && stop_loss_level)
                     if(!m_trade.PositionClose(m_position.Ticket())) // close a position by the specified m_symbol
                        Print(", ERROR: ","CTrade.PositionClose ",m_position.Ticket());
                 }
              }
  }
//+------------------------------------------------------------------+
//| Close Profitable positions                                       |
//+------------------------------------------------------------------+
void ClosePositionsProfitable(const ENUM_POSITION_TYPE pos_type,const double level,const double eThr)
  {
   /*
      Buying is done at the Ask price                 |  Selling is done at the Bid price
      ------------------------------------------------|----------------------------------
      TakeProfit        >= Bid                        |  TakeProfit        <= Ask
      StopLoss          <= Bid                        |  StopLoss          >= Ask
      TakeProfit - Bid  >= SYMBOL_TRADE_STOPS_LEVEL   |  Ask - TakeProfit  >= SYMBOL_TRADE_STOPS_LEVEL
      Bid - StopLoss    >= SYMBOL_TRADE_STOPS_LEVEL   |  StopLoss - Ask    >= SYMBOL_TRADE_STOPS_LEVEL
   */
   for(int i=PositionsTotal()-1; i>=0; i--) // returns the number of current positions
      if(m_position.SelectByIndex(i)) // selects the position by index for further access to its properties
         if(m_position.Symbol()==m_symbol.Name() && m_position.Magic()==InpMagic)
            if(m_position.PositionType()==pos_type && m_position.Profit()>=(eThr-m_position.Swap()-m_position.Commission()))
              {
               if(m_position.PositionType()==POSITION_TYPE_BUY)
                 {
                  bool take_profit_level=(m_position.TakeProfit()!=0.0 && m_position.TakeProfit()-m_position.PriceCurrent()>=level) || m_position.TakeProfit()==0.0;
                  bool stop_loss_level=(m_position.StopLoss()!=0.0 && m_position.PriceCurrent()-m_position.StopLoss()>=level) || m_position.StopLoss()==0.0;
                  if(take_profit_level && stop_loss_level)
                     if(!m_trade.PositionClose(m_position.Ticket())) // close a position by the specified m_symbol
                        Print(", ERROR: ","CTrade.PositionClose ",m_position.Ticket());
                 }
               if(m_position.PositionType()==POSITION_TYPE_SELL)
                 {
                  bool take_profit_level=(m_position.TakeProfit()!=0.0 && m_position.PriceCurrent()-m_position.TakeProfit()>=level) || m_position.TakeProfit()==0.0;
                  bool stop_loss_level=(m_position.StopLoss()!=0.0 && m_position.StopLoss()-m_position.PriceCurrent()>=level) || m_position.StopLoss()==0.0;
                  if(take_profit_level && stop_loss_level)
                     if(!m_trade.PositionClose(m_position.Ticket())) // close a position by the specified m_symbol
                        Print(", ERROR: ","CTrade.PositionClose ",m_position.Ticket());
                 }
              }
  }
//+------------------------------------------------------------------+
//| Close Dangerous positions                                        |
//+------------------------------------------------------------------+
void ClosePositionsDangers(const ENUM_POSITION_TYPE pos_type,const double level,const double eThr)
  {
   /*
      Buying is done at the Ask price                 |  Selling is done at the Bid price
      ------------------------------------------------|----------------------------------
      TakeProfit        >= Bid                        |  TakeProfit        <= Ask
      StopLoss          <= Bid                        |  StopLoss          >= Ask
      TakeProfit - Bid  >= SYMBOL_TRADE_STOPS_LEVEL   |  Ask - TakeProfit  >= SYMBOL_TRADE_STOPS_LEVEL
      Bid - StopLoss    >= SYMBOL_TRADE_STOPS_LEVEL   |  StopLoss - Ask    >= SYMBOL_TRADE_STOPS_LEVEL
   */
   for(int i=PositionsTotal()-1; i>=0; i--) // returns the number of current positions
      if(m_position.SelectByIndex(i)) // selects the position by index for further access to its properties
         if(m_position.Symbol()==m_symbol.Name() && m_position.Magic()==InpMagic)
            if(m_position.PositionType()==pos_type && m_position.Profit()<(-eThr-m_position.Swap()-m_position.Commission()))
              {
               if(m_position.PositionType()==POSITION_TYPE_BUY)
                 {
                  bool take_profit_level=(m_position.TakeProfit()!=0.0 && m_position.TakeProfit()-m_position.PriceCurrent()>=level) || m_position.TakeProfit()==0.0;
                  bool stop_loss_level=(m_position.StopLoss()!=0.0 && m_position.PriceCurrent()-m_position.StopLoss()>=level) || m_position.StopLoss()==0.0;
                  if(take_profit_level && stop_loss_level)
                     if(!m_trade.PositionClose(m_position.Ticket())) // close a position by the specified m_symbol
                        Print(", ERROR: ","CTrade.PositionClose ",m_position.Ticket());
                 }
               if(m_position.PositionType()==POSITION_TYPE_SELL)
                 {
                  bool take_profit_level=(m_position.TakeProfit()!=0.0 && m_position.PriceCurrent()-m_position.TakeProfit()>=level) || m_position.TakeProfit()==0.0;
                  bool stop_loss_level=(m_position.StopLoss()!=0.0 && m_position.StopLoss()-m_position.PriceCurrent()>=level) || m_position.StopLoss()==0.0;
                  if(take_profit_level && stop_loss_level)
                     if(!m_trade.PositionClose(m_position.Ticket())) // close a position by the specified m_symbol
                        Print(", ERROR: ","CTrade.PositionClose ",m_position.Ticket());
                 }
              }
  }
//+------------------------------------------------------------------+
//| Calculate all positions                                          |
//+------------------------------------------------------------------+
void CalculateAllPositions(int &count_buys,double &volume_buys,double &volume_biggest_buys,
                           int &count_sells,double &volume_sells,double &volume_biggest_sells,
                           int &count_hedge_buys, int &count_hedge_sells,
                           int &count_X_buys, int &count_X_sells)
  {
   count_buys  = 0;
   volume_buys   = 0.0;
   volume_biggest_buys  = 0.0;
   count_sells = 0;
   volume_sells  = 0.0;
   volume_biggest_sells = 0.0;
   for(int i=PositionsTotal()-1; i>=0; i--)
      if(m_position.SelectByIndex(i)) // selects the position by index for further access to its properties
         if(m_position.Symbol()==m_symbol.Name() && m_position.Magic()==InpMagic)
           {
            string PosComment = m_position.Comment();
            string PosType = StringSubstr(PosComment,0,10);
            long   PosIden = StringToInteger(StringSubstr(PosComment,11,-1));

            if(m_position.PositionType()==POSITION_TYPE_BUY)
              {
               if(PosType == "Normal****")
                  count_buys++;
               if(PosType == "Hedge_Buyy")
                  count_hedge_buys++;
               if(PosType == "XHedge_Buy")
                  count_X_buys++;
               volume_buys+=m_position.Volume();
               if(m_position.Volume()>volume_biggest_buys)
                  volume_biggest_buys=m_position.Volume();
               continue;
              }
            else
               if(m_position.PositionType()==POSITION_TYPE_SELL)
                 {
                  if(PosType == "Normal****")
                     count_sells++;
                  if(PosType == "Hedge_Sell")
                     count_hedge_sells++;
                  if(PosType == "XHedge_Sel")
                     count_X_sells++;
                  volume_sells+=m_position.Volume();
                  if(m_position.Volume()>volume_biggest_sells)
                     volume_biggest_sells=m_position.Volume();
                 }
           }
  }
//
//+------------------------------------------------------------------+
bool TimeControl(void)
  {
   if(!InpTimeControl)
      return(true);
   MqlDateTime STimeCurrent;
   datetime time_current=TimeCurrent();
   if(time_current==D'1970.01.01 00:00')
      return(false);
   TimeToStruct(time_current,STimeCurrent);
   if(InpStartHour<InpEndHour) // intraday time interval
     {
      /*
      Example:
      input uchar    InpStartHour      = 5;        // Start hour
      input uchar    InpEndHour        = 10;       // End hour
      0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15
      _  _  _  _  _  +  +  +  +  +  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  +  +  +  +  +  _  _  _  _  _  _
      */
      if(STimeCurrent.hour>=InpStartHour && STimeCurrent.hour<InpEndHour)
         return(true);
     }
   else
      if(InpStartHour>InpEndHour) // time interval with the transition in a day
        {
         /*
         Example:
         input uchar    InpStartHour      = 10;       // Start hour
         input uchar    InpEndHour        = 5;        // End hour
         0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15
         _  _  _  _  _  _  _  _  _  _  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  _  _  _  _  _  +  +  +  +  +  +
         */
         if(STimeCurrent.hour>=InpStartHour || STimeCurrent.hour<InpEndHour)
            return(true);
        }
      else
         return(false);
//---
   return(false);
  }
//+------------------------------------------------------------------+
//| Search trading signals                                           |
//+------------------------------------------------------------------+
bool SearchTradingSignals(void)
  {
//--- update fractals level
   double fractals_upper[],fractals_lower[];
   MqlRates rates[];

   ArraySetAsSeries(fractals_upper,true);
   ArraySetAsSeries(fractals_lower,true);
   ArraySetAsSeries(rates,true);

   int start_pos=0;

   if(!iGetArray(handle_iFractals,UPPER_LINE,start_pos,xcount,fractals_upper) ||
      !iGetArray(handle_iFractals,LOWER_LINE,start_pos,xcount,fractals_lower) ||
      CopyRates(m_symbol.Name(),Period(),start_pos,128,rates)!=128)
     {
      return(false);
     }
//Print("!!! Start");
//---------Last Nfr Fractals-----------------------
   ArrayInitialize(fractals_up,EMPTY_VALUE);
   ArrayInitialize(fractals_down,EMPTY_VALUE);
   ArrayInitialize(fractals_up_number,-1);
   ArrayInitialize(fractals_down_number,-1);

   int kku = 0;
   int kkd = 0;
   for(int i=2; i<xcount && (kkd<Nfr || kku<Nfr); i++)
     {
      if(fractals_upper[i]!=0.0 && fractals_upper[i]!=EMPTY_VALUE)
         if(kku<Nfr && fractals_up[kku]==EMPTY_VALUE)
           {
            fractals_up[kku]=fractals_upper[i];
            fractals_up_number[kku]=i;
            kku++;
           }
      if(fractals_lower[i]!=0.0 && fractals_lower[i]!=EMPTY_VALUE)
         if(kkd<Nfr && fractals_down[kkd]==EMPTY_VALUE)
           {
            fractals_down[kkd]=fractals_lower[i];
            fractals_down_number[kkd]=i;
            kkd++;
           }
     }

   if(fractals_up[0]==EMPTY_VALUE || fractals_down[0]==EMPTY_VALUE)
     {
      Print("!!! No Fractal!!!");
      return(false);
     }

   int size_need_position=ArraySize(SPosition);
//-----------
   static double prev_frac_up = -1;
   static bool fractal_ok_buy = false;
   static bool ch_frac_up = false;
//-----------
   static double prev_frac_down = 1e6;
   static bool fractal_ok_sell = false;
   static bool ch_frac_down = false;
//-----------
   static double fractals_upf = 0;
   static double fractals_downf = 0;

   if(MathAbs(fractals_up[0]-prev_frac_up)>=DelFrac*m_symbol.Point())
     {
      ch_frac_up = true;
     }
   prev_frac_up = fractals_up[0];
   FracsUpMA.Add(fractals_up[0]);
   if(FracsUpMA.Total()>FractFilt)
      FracsUpMA.Delete(0);
   fractals_upf = EKOMeanX(FracsUpMA);

   if(MathAbs(fractals_down[0]-prev_frac_down)>=DelFrac*m_symbol.Point())
     {
      ch_frac_down = true;
     }
   prev_frac_down = fractals_down[0];
   FracsDownMA.Add(fractals_down[0]);
   if(FracsDownMA.Total()>FractFilt)
      FracsDownMA.Delete(0);
   fractals_downf = EKOMeanX(FracsDownMA);

//Print("!!!Fu: ",fractals_up[0],", FracUpMA: ",fractals_upf);
//Print("!!!Fd: ",fractals_down[0],", FracDownMA: ",fractals_downf);
//Print("!!!-------------------");
//-----------
//-----------
//-----------
   double mean_price = (rates[1].close+rates[1].high+rates[1].low)/3;
//-------------------------------------------------|
//-------------------Statistics--------------------|
//-------------------------------------------------|
   static double SrPerc=0;
   if(m_stat.Calculate())
     {
      SRDiff += (m_stat.SharpeRatio()-SRref)<0?MathAbs(m_stat.SharpeRatio()-SRref):0;
      SrPerc = MathTanh((m_stat.SharpeRatio()-SRref)/SRref);
     }
   else
      Print(m_stat.GetLastErrorString());
//-----------
//-----------

//-----------
//-----------
   int count_buys=0;
   double volume_buys=0.0;
   double volume_biggest_buys=0.0;
   int count_sells=0;
   double volume_sells=0.0;
   double volume_biggest_sells=0.0;
   int    count_hedge_buys     = 0;
   int    count_hedge_sells    = 0;
   int    count_X_buys = 0;
   int    count_X_sells = 0;
   CalculateAllPositions(count_buys,volume_buys,volume_biggest_buys,
                         count_sells,volume_sells,volume_biggest_sells,
                         count_hedge_buys,count_hedge_sells,
                         count_X_buys,count_X_sells);
//-------------------------------------------------|
//-----------------------Buy-----------------------|
//-------------------------------------------------|
//---------------------
   datetime time_now = TimeCurrent();
   MqlDateTime date_now;
   TimeToStruct(time_now, date_now);
   static uint hournow = 1000;

   if(date_now.hour % 2 == 0 && date_now.hour != hournow && JustBuys)
     {
      OK_Buy = true;
      OK_Buy_Count = 0;
      hournow = date_now.hour;
     }

   OK_Buy_Count++;
   if(date_now.hour % 2 == 1 ||
      OK_Buy_Count>OKAcceptCount)
      OK_Buy_Count = OKAcceptCount+1;

   if((OK_Buy && OK_Buy_Count<OKAcceptCount) &&/*
      ((m_symbol.Ask()/fractals_upf)>(1+fraka/100)) &&
      (ch_frac_down) &&*/
      ((InpReverse?count_sells:count_buys)<((shadowEn?2:1)*ttd)))
     {
      //Print("!!!Frac_Buy: ",m_symbol.Ask(),">",fractals_down[0]);
      if(InpReverse)
         Print("!!!Count Sells: ",count_sells);
      else
         Print("!!!Count Buys: ",count_buys);
      Print("!!!Sharpe Ratio: ",m_stat.SharpeRatio());
      Print("!!!LR Correlation: ",m_stat.LRCorrelation());
      Print("!!!Z-Score: ",m_stat.ZScore());
      Print("!!!-----------------------------");
      CountTrades++;
      if(CountTrades>InitTradeCount)
         CountTrades = InitTradeCount;
      if(shadowEn)
        {
         ArrayResize(SPosition,size_need_position+2);
         SPosition[size_need_position].pos_type=TradeX;
         SPosition[size_need_position].comment = "Normal****";
         if(BrakeEn)
           {
            SPosition[size_need_position].lot_coefficient *= BrVolume;
           }

         if(CountTrades<InitTradeCount && InpLotOrRisk==risk)
            SPosition[size_need_position].volume = m_symbol.LotsMin();
         if(revShadow)
            SPosition[size_need_position+1].pos_type=rTradeX;
         else
            SPosition[size_need_position+1].pos_type=TradeX;
         SPosition[size_need_position+1].comment = "Normal****";
         SPosition[size_need_position+1].stopLossPercent = slPercent;
         SPosition[size_need_position+1].takeProfitPercent = tpPercent;
         if(BrakeEn)
           {
            SPosition[size_need_position+1].lot_coefficient *= BrVolume;
           }

         if(CountTrades<InitTradeCount && InpLotOrRisk==risk)
            SPosition[size_need_position+1].volume = m_symbol.LotsMin();
        }
      else
        {
         ArrayResize(SPosition,size_need_position+1);
         SPosition[size_need_position].pos_type=TradeX;
         SPosition[size_need_position].comment = "Normal****";
         if(BrakeEn)
           {
            SPosition[size_need_position].lot_coefficient *= BrVolume;
           }

         if(CountTrades<InitTradeCount && InpLotOrRisk==risk)
            SPosition[size_need_position].volume = m_symbol.LotsMin();
        }
      //---------------------
      ch_frac_down = false;
      countBuys++;
      OK_Buy = false;
     }

   CalculateAllPositions(count_buys,volume_buys,volume_biggest_buys,
                         count_sells,volume_sells,volume_biggest_sells,
                         count_hedge_buys,count_hedge_sells,
                         count_X_buys,count_X_sells);
   if((InpReverse?count_sells:count_buys)>0 && enExit && JustBuys)
     {
      double level;
      double X1 = m_account.Balance();
      if(FreezeStopsLevels(level))
        {
         ClosePositionsProfitable(TradeX,level,eThreshold);
        }
      if((m_account.Balance()-X1)>0)
        {
         countClose++;
        }
     }
//-------------------------------------------------|
//-------------------End of Buy--------------------|
//-------------------------------------------------|
   CalculateAllPositions(count_buys,volume_buys,volume_biggest_buys,
                         count_sells,volume_sells,volume_biggest_sells,
                         count_hedge_buys,count_hedge_sells,
                         count_X_buys,count_X_sells);
//-------------------------------------------------|
//---------------------Sell------------------------|
//-------------------------------------------------|
   if(date_now.hour % 2 == 1 && date_now.hour != hournow && JustSells)
     {
      OK_Sell = true;
      OK_Sell_Count = 0;
      hournow = date_now.hour;
     }

   OK_Sell_Count++;
   if(date_now.hour % 2 == 0 ||
      OK_Sell_Count>OKAcceptCount)
      OK_Sell_Count = OKAcceptCount+1;

   if((OK_Sell && OK_Sell_Count<OKAcceptCount) &&/*
      ((m_symbol.Bid()/fractals_downf)<(1-fraka/100)) &&
      (ch_frac_up) &&*/
      ((InpReverse?count_buys:count_sells)<((shadowEn?2:1)*ttd)))
     {
      //Print("!!!Frac_Sell: ",m_symbol.Bid(),"<",fractals_up[0]);
      if(InpReverse)
         Print("!!!Count Buys: ",count_buys);
      else
         Print("!!!Count Sells: ",count_sells);

      Print("!!!Sharpe Ratio: ",m_stat.SharpeRatio());
      Print("!!!LR Correlation: ",m_stat.LRCorrelation());
      Print("!!!Z-Score: ",m_stat.ZScore());
      Print("!!!-----------------------------");
      CountTrades++;
      if(CountTrades>InitTradeCount)
         CountTrades = InitTradeCount;
      if(shadowEn)
        {
         ArrayResize(SPosition,size_need_position+2);
         SPosition[size_need_position].pos_type=rTradeX;
         SPosition[size_need_position].comment = "Normal****";
         if(BrakeEn)
           {
            SPosition[size_need_position].lot_coefficient *= BrVolume;
           }

         if(CountTrades<InitTradeCount && InpLotOrRisk==risk)
            SPosition[size_need_position].volume = m_symbol.LotsMin();

         if(revShadow)
            SPosition[size_need_position+1].pos_type=TradeX;
         else
            SPosition[size_need_position+1].pos_type=rTradeX;
         SPosition[size_need_position+1].comment = "Normal****";
         SPosition[size_need_position+1].stopLossPercent = slPercent;
         SPosition[size_need_position+1].takeProfitPercent = tpPercent;
         if(BrakeEn)
           {
            SPosition[size_need_position+1].lot_coefficient *= BrVolume;
           }

         if(CountTrades<InitTradeCount && InpLotOrRisk==risk)
            SPosition[size_need_position+1].volume = m_symbol.LotsMin();
        }
      else
        {
         ArrayResize(SPosition,size_need_position+1);
         SPosition[size_need_position].pos_type=rTradeX;
         SPosition[size_need_position].comment = "Normal****";
         if(BrakeEn)
           {
            SPosition[size_need_position].lot_coefficient *= BrVolume;
           }

         if(CountTrades<InitTradeCount && InpLotOrRisk==risk)
            SPosition[size_need_position].volume = m_symbol.LotsMin();
        }

      ch_frac_up = false;
      countSells++;
      OK_Sell = false;
     }

   CalculateAllPositions(count_buys,volume_buys,volume_biggest_buys,
                         count_sells,volume_sells,volume_biggest_sells,
                         count_hedge_buys,count_hedge_sells,
                         count_X_buys,count_X_sells);
   if((InpReverse?count_buys:count_sells)>0 && enExit && JustSells)
     {
      double level;
      double X1 = m_account.Balance();
      if(FreezeStopsLevels(level))
        {
         ClosePositionsProfitable(rTradeX,level,eThreshold);
        }
      if((m_account.Balance()-X1)>0)
        {
         countClose++;
        }
     }
//-------------------------------------------------|
//--------------end of Sell------------------------|
//-------------------------------------------------|
   return(true);
  }

//+------------------------------------------------------------------+
//|----------------------Array Normalization-------------------------|
//+------------------------------------------------------------------+
void NormalizeArrays(double &a[])
  {
   double d1=0;
   double d2=1.0;
   double x_min=0.0;
   double x_max=100;
   for(int i=0; i<ArraySize(a); i++)
     {
      a[i]=(((a[i]-x_min)*(d2-d1))/(x_max-x_min))+d1;
     }
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
double EKOMean(CArrayDouble &Cinp)
  {
   int size = Cinp.Total();
   double sum = 0;
   for(int i=0; i<size; ++i)
     {
      sum+=Cinp.At(i);
     }
   return sum/size;
  }

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
double EKOMeanX(CArrayDouble &Cinp)
  {
   int size = Cinp.Total();
   double sum = 0;
   for(int i=0; i<size; ++i)
     {
      sum+=(i+1)*(i+1)*Cinp.At(i);
     }
   return sum/(size*(size+1)*(2*size+1)/6);
  }
//+------------------------------------------------------------------+
//|----------------------Expiration----------------------------------|
//+------------------------------------------------------------------+
bool EKOExpire(const ENUM_POSITION_TYPE pos_type,const double level,const double eThr,const long ExpMinute)
  {
   if(ExpireEn)
     {
      for(int i=PositionsTotal()-1; i>=0; i--) // returns the number of current positions
         if(m_position.SelectByIndex(i)) // selects the position by index for further access to its properties
            if(m_position.Symbol()==m_symbol.Name() && m_position.Magic()==InpMagic)
               if(CheckTime(m_position.Time(),ExpMinute))
                  if(m_position.PositionType()==pos_type &&
                     m_position.Profit()>=(eThr-m_position.Swap()-m_position.Commission()))
                    {
                     Print("!!!Position Expired; Swap: ",m_position.Swap(),"; Commision: ",m_position.Commission());
                     if(m_position.PositionType()==POSITION_TYPE_BUY)
                       {
                        if(!m_trade.PositionClose(m_position.Ticket())) // close a position by the specified m_symbol
                           Print(", ERROR: ","CTrade.PositionClose ",m_position.Ticket());
                       }
                     if(m_position.PositionType()==POSITION_TYPE_SELL)
                       {
                        if(!m_trade.PositionClose(m_position.Ticket())) // close a position by the specified m_symbol
                           Print(", ERROR: ","CTrade.PositionClose ",m_position.Ticket());
                       }
                    }
     }
   return true;
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
bool EKOExpireTerm(const ENUM_POSITION_TYPE pos_type,const double level,const double TermThr,const long ExpMinute)
  {
   if(ExTermEn)
     {
      for(int i=PositionsTotal()-1; i>=0; i--) // returns the number of current positions
         if(m_position.SelectByIndex(i)) // selects the position by index for further access to its properties
            if(m_position.Symbol()==m_symbol.Name() && m_position.Magic()==InpMagic)
              {
               string PosComment = m_position.Comment();
               string PosType = StringSubstr(PosComment,0,10);
               long   PosIden = StringToInteger(StringSubstr(PosComment,11,-1));

               if(PosType=="Normal****" && CheckTime(m_position.Time(),2*ExpMinute))
                  if(m_position.PositionType()==pos_type &&
                     MathAbs(m_account.Balance()-m_account.Equity())<=TermThr)
                    {
                     Print("!!!Position Terminated at: ",MathAbs(m_account.Balance()-m_account.Equity()),
                        "; Swap: ",m_position.Swap(),"; Commision: ",m_position.Commission());
                     if(m_position.PositionType()==POSITION_TYPE_BUY)
                       {
                        if(!m_trade.PositionClose(m_position.Ticket())) // close a position by the specified m_symbol
                           Print(", ERROR: ","CTrade.PositionClose ",m_position.Ticket());
                       }
                     if(m_position.PositionType()==POSITION_TYPE_SELL)
                       {
                        if(!m_trade.PositionClose(m_position.Ticket())) // close a position by the specified m_symbol
                           Print(", ERROR: ","CTrade.PositionClose ",m_position.Ticket());
                       }
                    }
              }
     }
   return true;
  }
//+------------------------------------------------------------------+
bool EKOBraked()
  {
   for(int i=PositionsTotal()-1; i>=0; i--) // returns the number of current positions
      if(m_position.SelectByIndex(i)) // selects the position by index for further access to its properties
         if(m_position.Symbol()==m_symbol.Name() && m_position.Magic()==InpMagic)
           {
            string PosComment = m_position.Comment();
            string PosType = StringSubstr(PosComment,0,10);
            long   PosIden = StringToInteger(StringSubstr(PosComment,11,-1));

            if(PosType=="Normal****")
              {
               //Print("!!!Position Braked");
               if(!m_trade.PositionClose(m_position.Ticket())) // close a position by the specified m_symbol
                  Print(", ERROR: ","CTrade.PositionClose ",m_position.Ticket());
              }
           }
   return true;
  }

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
bool EKOBrakedX(double eth)
  {
   for(int i=PositionsTotal()-1; i>=0; i--) // returns the number of current positions
      if(m_position.SelectByIndex(i)) // selects the position by index for further access to its properties
         if(m_position.Symbol()==m_symbol.Name() && m_position.Magic()==InpMagic)
           {
            string PosComment = m_position.Comment();
            string PosType = StringSubstr(PosComment,0,10);
            long   PosIden = StringToInteger(StringSubstr(PosComment,11,-1));

            if(PosType=="Normal****" &&
               m_position.Profit()<=(-eth-m_position.Swap()-m_position.Commission()))
              {
               //Print("!!!Position Braked");
               if(!m_trade.PositionClose(m_position.Ticket())) // close a position by the specified m_symbol
                  Print(", ERROR: ","CTrade.PositionClose ",m_position.Ticket());
              }
           }
   return true;
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
bool CheckTime(datetime inTime, long compTimeMinute)
  {
   long time_now = (long) TimeCurrent();
//return ((time_now-(long)inTime)/60)>compTimeMinute?true:false;
   if(((time_now-(long)inTime)/60)>compTimeMinute)
     {
      sleep_time++;
      return true;
     }
   else
      return false;
  }
//+------------------------------------------------------------------+
//|----------------------Hedging-------------------------------------|
//+------------------------------------------------------------------+
bool EKOHedge(void)
  {
   if(HedgeEn)
     {
      //Print("!!!Enter Hedge");
      int count_buys=0;
      double volume_buys=0.0;
      double volume_biggest_buys=0.0;
      int count_sells=0;
      double volume_sells=0.0;
      double volume_biggest_sells=0.0;
      int    count_hedge_buys     = 0;
      int    count_hedge_sells    = 0;
      int    count_X_buys     = 0;
      int    count_X_sells    = 0;

      CalculateAllPositions(count_buys,volume_buys,volume_biggest_buys,
                            count_sells,volume_sells,volume_biggest_sells,
                            count_hedge_buys,count_hedge_sells,
                            count_X_buys,
                            count_X_sells);

      HedgedPositionsID.Sort();
      for(int i=PositionsTotal()-1; i>=0; i--)
         if(m_position.SelectByIndex(i)) // selects the position by index for further access to its properties
            if(m_position.Symbol()==m_symbol.Name() && m_position.Magic()== InpMagic)
              {
               if(m_position.PositionType()==POSITION_TYPE_BUY)
                 {
                  string PosComment = m_position.Comment();
                  string PosType = StringSubstr(PosComment,0,10);
                  long   PosIden = StringToInteger(StringSubstr(PosComment,11,-1));

                  if(PosType=="Normal****" &&
                     (HedgedPositionsID.Search(m_position.Identifier())==-1))
                    {
                     if(((m_position.PriceCurrent()-m_position.PriceOpen())<=-Hedge_DrawdownOpen))
                       {
                        int xsize_need_position=ArraySize(SPosition);
                        ArrayResize(SPosition,xsize_need_position+1);
                        SPosition[xsize_need_position].pos_type= POSITION_TYPE_SELL;
                        SPosition[xsize_need_position].comment = "Hedge_Sell,"+IntegerToString(m_position.Identifier());
                        SPosition[xsize_need_position].volume = Htype? m_position.Volume(): m_position.Volume()*Hedge_Coef*((count_X_buys+1));
                        SPosition[xsize_need_position].stopLossPercent = 1;
                        SPosition[xsize_need_position].takeProfitPercent = 1;
                        SPosition[xsize_need_position].Identifier = m_position.Identifier();
                        HedgedPositionsID.Add(SPosition[xsize_need_position].Identifier);

                        if(ReHedgeEn)
                          {
                           xsize_need_position=ArraySize(SPosition);
                           ArrayResize(SPosition,xsize_need_position+1);
                           SPosition[xsize_need_position].pos_type= POSITION_TYPE_BUY;
                           SPosition[xsize_need_position].comment = "XHedge_Buy,"+IntegerToString(m_position.Identifier());
                           SPosition[xsize_need_position].volume = m_position.Volume();
                           SPosition[xsize_need_position].stopLossPercent = 1;
                           SPosition[xsize_need_position].takeProfitPercent = 1;
                           SPosition[xsize_need_position].Identifier = m_position.Identifier();
                           countXBuys++;
                          }
                       }
                    }
                  if(PosType=="Hedge_Buyy")
                    {
                     if(m_position.SelectByTicket(PosIden))
                       {
                        double Xprof = (m_position.PriceCurrent()-m_position.PriceOpen());
                        m_position.SelectByIndex(i);
                        if((Xprof)<-Hedge_DrawdownOpen ||
                           ((m_position.PriceCurrent()-m_position.PriceOpen())>Hedge_Min_Profit && m_position.Profit()>0))
                          {
                           if(m_trade.PositionClose(m_position.Ticket()))
                              HedgedPositionsID.Delete(HedgedPositionsID.Search(PosIden));
                          }
                       }
                     else
                       {
                        m_position.SelectByIndex(i);
                        if(m_trade.PositionClose(m_position.Ticket()))
                           HedgedPositionsID.Delete(HedgedPositionsID.Search(PosIden));
                       }
                    }
                  if(PosType=="XHedge_Buy")
                    {
                     if(!m_position.SelectByTicket(PosIden))
                       {
                        m_position.SelectByIndex(i);
                        m_trade.PositionClose(m_position.Ticket());
                       }
                    }
                 }//Buy
               else
                  if(m_position.PositionType()==POSITION_TYPE_SELL)
                    {
                     string PosComment = m_position.Comment();
                     string PosType = StringSubstr(PosComment,0,10);
                     long   PosIden = StringToInteger(StringSubstr(PosComment,11,-1));

                     if(PosType=="Normal****" &&
                        (HedgedPositionsID.Search(m_position.Identifier())==-1))
                       {
                        if(((m_position.PriceCurrent()-m_position.PriceOpen())>=Hedge_DrawdownOpen))
                          {
                           int xsize_need_position=ArraySize(SPosition);
                           ArrayResize(SPosition,xsize_need_position+1);
                           SPosition[xsize_need_position].pos_type= POSITION_TYPE_BUY;
                           SPosition[xsize_need_position].comment = "Hedge_Buyy,"+IntegerToString(m_position.Identifier());
                           SPosition[xsize_need_position].volume = Htype? m_position.Volume(): m_position.Volume()*Hedge_Coef*((count_X_sells+1));
                           SPosition[xsize_need_position].stopLossPercent = 1;
                           SPosition[xsize_need_position].takeProfitPercent = 1;
                           SPosition[xsize_need_position].Identifier = m_position.Identifier();
                           HedgedPositionsID.Add(SPosition[xsize_need_position].Identifier);

                           if(ReHedgeEn)
                             {
                              xsize_need_position=ArraySize(SPosition);
                              ArrayResize(SPosition,xsize_need_position+1);
                              SPosition[xsize_need_position].pos_type= POSITION_TYPE_SELL;
                              SPosition[xsize_need_position].comment = "XHedge_Sel,"+IntegerToString(m_position.Identifier());
                              SPosition[xsize_need_position].volume = m_position.Volume();
                              SPosition[xsize_need_position].stopLossPercent = 1;
                              SPosition[xsize_need_position].takeProfitPercent = 1;
                              SPosition[xsize_need_position].Identifier = m_position.Identifier();
                              countXSells++;
                             }
                          }
                       }
                     if(PosType=="Hedge_Sell")
                       {
                        if(m_position.SelectByTicket(PosIden))
                          {
                           double Xprof = (m_position.PriceCurrent()-m_position.PriceOpen());
                           m_position.SelectByIndex(i);
                           if((Xprof)>Hedge_DrawdownOpen ||
                              ((m_position.PriceCurrent()-m_position.PriceOpen())<-Hedge_Min_Profit && m_position.Profit()>0))
                             {
                              if(m_trade.PositionClose(m_position.Ticket()))
                                 HedgedPositionsID.Delete(HedgedPositionsID.Search(PosIden));
                             }
                          }
                        else
                          {
                           m_position.SelectByIndex(i);
                           if(m_trade.PositionClose(m_position.Ticket()))
                              HedgedPositionsID.Delete(HedgedPositionsID.Search(PosIden));
                          }
                       }
                     if(PosType=="XHedge_Sel")
                       {
                        if(!m_position.SelectByTicket(PosIden))
                          {
                           m_position.SelectByIndex(i);
                           m_trade.PositionClose(m_position.Ticket());
                          }
                       }
                    }//Sell
              }//Magic
     }

   return(true);
  }
//-------------------------------------------------|
//--------------------Hedging End------------------|
//-------------------------------------------------|
double ThrX(double inpValue, double alp, double Tr)
  {
   double x = alp*inpValue;
   x = x>Tr?Tr:x;
   x = x<-Tr?-Tr:x;
   return x;
  }

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
double EKOSign(double x, double tr)
  {
   if(x>=tr)
      return 1;
   if(MathAbs(x)<tr)
      return 0;
   if(x<=-tr)
      return -1;
   return 0;
  }
//--------------------------------------------------------------------
double OnTester()
  {
   double  param = 0.0;

//  Balance max + min Drawdown + Trades Number:
   static double  eP = TesterStatistics(STAT_EXPECTED_PAYOFF);
   static double  bl = TesterStatistics(STAT_PROFIT);
   static double  bS = TesterStatistics(STAT_GROSS_LOSS);
   static double  pF = TesterStatistics(STAT_PROFIT_FACTOR);
   static double  dED = TesterStatistics(STAT_EQUITYDD_PERCENT);
   static double  dBD = TesterStatistics(STAT_BALANCEDD_PERCENT);
   static double  rF = TesterStatistics(STAT_RECOVERY_FACTOR);
   static double  sR = TesterStatistics(STAT_SHARPE_RATIO);
   static double  tP = TesterStatistics(STAT_PROFIT_TRADES);
   static double  tL = TesterStatistics(STAT_LOSS_TRADES);
   static double  tT  = TesterStatistics(STAT_TRADES);
   static double  mL  = TesterStatistics(STAT_MIN_MARGINLEVEL);
   int mm=1;
   double pFlim = 2;
   mL = mL>0?mL:0;

   double  rED = 1;
   if(dED > 0.01)
     {
      rED = 1.0 / dED;
     }
   else
     {
      rED = 0;
     }

   double  rBD = 1;
   if(dBD > 0.01)
     {
      rBD = 1.0 / dBD;
     }
   else
     {
      rBD = 0;
     }

   double rBl = bl>0?bl:0;
   if(rBl<BlLim)
     {
      bl = rBl*MathExp(-(BlLim-rBl)/5);
     }

   if(rF<0 || rF>32)
     {
      rF = 0;
     }

   if(pF<0)
     {
      pF = 0;
     }

   if(eP<0)
     {
      eP = 0;
     }

   double rSr = sR>0?sR:0;
   if(sR<=SharpeLow || sR>(SharpeLow+SharpeHighPlus))
     {
      sR = 0;
     }

   if(tP<5)
      tP = 0;

   double mntcrlo = optpr();

   if(MWE)
     {
      param = MathPow(wk*co_money,BlFac)*tP*tT*MathPow(sR,SRFac)*SRFac*DWFactor*(SharpeLow+0.1)*BlLim*mL*MathExp(SRref)*
              MathPow(eP*rF,EpRfFac)*MathSqrt((tT<(XtT/mm)?0:XtT*tT)*(countClose+1))*wd_money*(co_money>BlLim?MathExp(3*sR):1)*
              MathSqrt((rBD*rED+1))*MathSqrt(TermExtract)/

              MathSqrt(SharpeHighPlus*InpTakeProfit*InpStopLoss*MathExp(MathSqrt(InpVolumeLotOrRisk))*MathExp(ttd)*(BrakeEn?(MathExp(brakeCountReal)+1):1)*
                       (JustBuys?1:MathExp(MathAbs(countBuys-countSells)/MathLog(bl/2+2)))*ExpireMinute*(tL+1)*InpVolumeLotOrRisk*InpDeviation*
                       (tT<XtT?MathExp((XtT-tT)/MathLog(bl/2+2)):1/MathSqrt(tT*tT-XtT*XtT+1))*ExpireExtract*(Positivised+1)*gotoActiveCount*SRDiff*
                       MathPow(under_MaxEq*under_Linear+1,DWFactor)*MathExp(BrakePerT)*(MathAbs(BrakeThreshold)+1)*MathSqrt(TermExtract)+1)+

              (starterEn?MathLog(countSells*countBuys*rBl*tP*tT*rSr/(InpTakeProfit*InpStopLoss*MathExp(InpVolumeLotOrRisk/2)*
                                 MathExp(ttd+1)*(tL+1)*(Positivised+1)+1)+1):0+

               (starterEn?rBl*MathSqrt(countBuys*countSells)*(mntcrlo>0?mntcrlo:0)/
                InpVolumeLotOrRisk:0));

      param = MCOptEn?rBl*MathSqrt(tT)*(mntcrlo>0?mntcrlo:0)/
              MathSqrt(InpVolumeLotOrRisk*InpTakeProfit*InpStopLoss+1)/
              (JustBuys?1:MathExp(0.5*MathAbs(countBuys-countSells)/MathSqrt(tT+1))):param;

     }
   else
     {

      param = (pF<pFlim?0:MathSqrt(pF))*MathPow(bl,BlFac)*tP*tT*MathPow(sR,SRFac)*SRFac*DWFactor*(SharpeLow+0.1)*BlLim*mL*MathExp(SRref)*
              MathPow(eP*rF,EpRfFac)*MathSqrt((tT<(XtT/mm)?0:XtT*tT)*(countClose+1))*(bl>BlLim?MathExp(3*sR):1)*
              MathSqrt((rBD*rED+1))*MathSqrt(TermExtract)/

              MathSqrt(SharpeHighPlus*InpTakeProfit*InpStopLoss*MathExp(MathSqrt(InpVolumeLotOrRisk))*MathExp(ttd)*(BrakeEn?(MathExp(brakeCountReal)+1):1)*
                       (JustBuys?1:MathExp(MathAbs(countBuys-countSells)/MathLog(bl/2+2)))*ExpireMinute*(tL+1)*InpVolumeLotOrRisk*InpDeviation*
                       (tT<XtT?MathExp((XtT-tT)/MathLog(bl/2+2)):1/MathSqrt(tT*tT-XtT*XtT+1))*ExpireExtract*(Positivised+1)*gotoActiveCount*SRDiff*
                       MathPow(under_MaxEq*under_Linear+1,DWFactor)*MathExp(BrakePerT)*(MathAbs(BrakeThreshold)+1)*MathSqrt(TermExtract)+1)+

              (starterEn?((pF<2?0:MathSqrt(pF))*MathLog(countSells*countBuys*rBl*tP*((tT<(XtT/mm) || rBl<tT)?0:tT)*rSr/

                          (InpTakeProfit*InpStopLoss*MathExp(InpVolumeLotOrRisk/2)*
                           MathExp(ttd+1)*(tL+1)*(Positivised+1)+1)+1)+

                          rBl*MathSqrt((tT<(XtT/mm) || rBl<tT)?0:rSr*tT)*(mntcrlo>0?mntcrlo:0)*(bl>BlLim?MathExp(sR):1)*(pF<2?0:MathSqrt(pF))/
                          MathSqrt(InpVolumeLotOrRisk*InpTakeProfit*InpStopLoss+1)/
                          (JustBuys?1:MathExp(0.5*MathAbs(countBuys-countSells)/MathSqrt(tT+1)))):0);
      param = param/(EqExValue+1);

      param = MCOptEn?rBl*MathSqrt(tT*rSr)*(mntcrlo>0?mntcrlo:0)*(pF<2?0:MathSqrt(pF))*(eP<1.1?0:MathSqrt(eP))/
              MathSqrt(InpVolumeLotOrRisk*InpTakeProfit*InpStopLoss+1)/
              (JustBuys?1:MathExp(0.5*MathAbs(countBuys-countSells)/MathSqrt(tT+1))):param;

      param = MaxBal?((rBl<tT/*2.5*init_Eq*/?0:rBl)*tT*(pF<2?MathSqrt(pF):pF)*(eP<1.1?MathSqrt(eP):eP)*(rF<2?MathSqrt(rF):rF)/(tL+1)):param;
      //tT=tT<100?0:tT;
      //param = NNtrainer?(pF>2.1?MathSqrt(pF*rBl*rSr):0)*(((rBl>tT)?rBl:0)*MathLog(tT*tP+1)+

      //        tT*tP*MathLog(((rBl>tT)?rBl:0)+1))/
      //        ((tL+1)*InpTakeProfit*InpStopLoss*MathExp(InpVolumeLotOrRisk/sqrt(rBl*rSr))*
      //         MathSqrt(((rBl>tT)?rBl:0)+tT*tP+1)):param;

      //param = NNtrainer?(pF>2?(rSr>0.1?(tT>XtT?MathSqrt(tT*rBl)*pF*rSr*InpVolumeLotOrRisk/(tL+1):0):0):0):param;MathSqrt

      param = NNtrainer?((pF<2?0:MathSqrt(pF))*(rBl<tT?0:rBl)*(tT<(XtT)?0:tT*tT)/(EqExValue*tL+1)):param;

     }

//param = optpr()*rBl*eP*rF*rSr*MathSqrt((tT<(XtT/mm) || rBl<tT)?0:XtT*tT)/
//        (BrakePerT*InpVolumeLotOrRisk*(BrakeEn?MathLog(brakeCount+1):1)+1);

   return(param);
  }
//+------------------------------------------------------------------+
