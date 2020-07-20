//
// Created by km on 08/06/2020.
//
#include <ql/cashflows/fixedratecoupon.hpp>
#include <ql/cashflows/iborcoupon.hpp>
#include <ql/currencies/europe.hpp>
#include <ql/indexes/ibor/euribor.hpp>
#include <ql/instruments/bonds/fixedratebond.hpp>
#include <ql/instruments/claim.hpp>
#include <ql/instruments/creditdefaultswap.hpp>
#include <ql/experimental/credit/riskybond.hpp>
#include <ql/currency.hpp>

#include <ql/math/interpolations/backwardflatinterpolation.hpp>
#include <ql/pricingengines/bond/bondfunctions.hpp>
#include <ql/pricingengines/bond/discountingbondengine.hpp>
#include <ql/pricingengines/credit/isdacdsengine.hpp>
#include <ql/pricingengines/credit/midpointcdsengine.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/credit/defaultprobabilityhelpers.hpp>
#include <ql/termstructures/credit/flathazardrate.hpp>
#include <ql/termstructures/credit/piecewisedefaultcurve.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/yield/piecewiseyieldcurve.hpp>
#include <ql/termstructures/yield/ratehelpers.hpp>
#include <ql/termstructures/yield/zerospreadedtermstructure.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/calendars/weekendsonly.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/time/daycounters/actualactual.hpp>
#include <ql/time/daycounters/thirty360.hpp>
#include <iomanip>
#include <iostream>

using namespace std;
using namespace QuantLib;

void example_cds()
{
    Calendar calendar = TARGET() ;
    Date tradeDate(12, May, 2020) ;
    Settings::instance().evaluationDate() = tradeDate;

    // flat riskless discounting
    // dummy curve
    ext::shared_ptr<Quote> flatRate(new SimpleQuote(0.01));
    Handle<YieldTermStructure> tsCurve(
        ext::make_shared<FlatForward>(tradeDate,
                                      Handle<Quote>(flatRate),
                                      Actual365Fixed()));


    // market
    Real recovery_rate = 0.5;
    Real quoted_spreads[] = { 0.0150, 0.0150, 0.0150, 0.0150 };
    vector<Period> tenors = {3 * Months, 6 * Months, 1 * Years, 2 * Years};

    tenors.push_back(3 * Months);
    tenors.push_back(6 * Months);
    tenors.push_back(1 * Years);
    tenors.push_back(2 * Years);
    vector<Date> maturities;
    for (Size i = 0; i < 4; i++) {
        maturities.push_back(
            calendar.adjust(tradeDate + tenors[i], Following));
//        std::cout<<"Maturity: "<<std::endl;
//        std::cout<<maturities[i]<<endl ;
    }

    std::vector<ext::shared_ptr<DefaultProbabilityHelper> > instruments;
    for (Size i = 0; i < 4; i++) {
        instruments.push_back(ext::shared_ptr<DefaultProbabilityHelper>(
            new SpreadCdsHelper(Handle<Quote>(ext::shared_ptr<Quote>(
                new SimpleQuote(quoted_spreads[i]))),
                                tenors[i], 0, calendar, Quarterly, Following,
                                DateGeneration::TwentiethIMM, Actual365Fixed(),
                                recovery_rate, tsCurve)));

    }

    // Bootstrap hazard rates
    ext::shared_ptr<PiecewiseDefaultCurve<HazardRate, BackwardFlat> >
        hazardRateStructure(new PiecewiseDefaultCurve<HazardRate, BackwardFlat>(tradeDate, instruments, Actual365Fixed()));

    Size UPFRONT_SETTLEMENT_DAYS(0);
    Handle<DefaultProbabilityTermStructure> probability(hazardRateStructure);
    ext::shared_ptr<PricingEngine> engine(
        new IsdaCdsEngine(probability, recovery_rate, tsCurve));

    Date upfrontDate = calendar.advance(tradeDate, UPFRONT_SETTLEMENT_DAYS * Period(Daily)) ;
    Real notional = 100000000. ;
    Date maturity_date(20, June, 2021);

    Schedule sched = MakeSchedule()
                     .from(tradeDate)
                     .to(maturity_date)
                     .withFrequency(Quarterly)
                     .withCalendar(calendar)
                     .withTerminationDateConvention(Unadjusted)
                     .withRule(QuantLib::DateGeneration::TwentiethIMM) ;


    ext::shared_ptr<Claim> claim_(new FaceValueClaim());

    CreditDefaultSwap cds_(Protection::Seller,
                           notional,
                           0.,
                           0.01,
                           sched,
                           Following,
                           Actual360(),
                           true,
                           true,
                           tradeDate,
                           upfrontDate,
                           claim_,
                           Actual360(true)
                           ) ;

    cds_.setPricingEngine(engine) ;
    cds_.fairUpfront() ;
}

std::vector<ext::shared_ptr<RateHelper>>
get_isda_curve_nodes(const std::vector<Real>& dep_rates,
                       const std::vector<Real>& swap_rates,
                       const Calendar& deposit_calendar,
                       const BusinessDayConvention& deposit_conv,
                       const DayCounter& deposit_day_counter,
                       bool end_of_month,
                       const Calendar& swCalendar,
                       const Frequency& swFixedLegFrequency,
                       const BusinessDayConvention& swFixedLegConv,
                       const DayCounter& swFixLegDayCounter,
                       ext::shared_ptr<IborIndex> swFltLegIndex)
{
    Integer isda_settlement_days = 3 ;

    std::vector<ext::shared_ptr<RateHelper>> depoAndSwaps ;
    depoAndSwaps.reserve(dep_rates.size() + swap_rates.size());

    std::vector<Period> dep_rates_tenors = {1 * Months, 3 * Months, 6 * Months, 12 * Months};
    std::vector<Period> swap_rates_tenors = {
        2 * Years,
        3 * Years,
        5 * Years,
        10 * Years,
        15 * Years,
    };

    for(Size ii = 0 ; ii < dep_rates.size() ; ii++) {
        ext::shared_ptr<Quote> dep_rate(new SimpleQuote(dep_rates[ii]));
        ext::shared_ptr<RateHelper> depRateHelper(new DepositRateHelper(
            Handle<Quote>(dep_rate), dep_rates_tenors[ii], isda_settlement_days, deposit_calendar,
            deposit_conv, end_of_month, deposit_day_counter));

        depoAndSwaps.push_back(depRateHelper);
    };

    for(Size ii = 0 ; ii < swap_rates.size() ; ii++) {
        ext::shared_ptr<Quote> swap_rate(new SimpleQuote(swap_rates[ii]));
        ext::shared_ptr<RateHelper> swRateHelper(new SwapRateHelper(
            Handle<Quote>(swap_rate),
            swap_rates_tenors[ii],
            swCalendar,
            swFixedLegFrequency,
            swFixedLegConv,
            swFixLegDayCounter,
            swFltLegIndex,
            Handle<Quote>(),
            Period(1 * Days)
            ));

        depoAndSwaps.push_back(swRateHelper);
    };
    return depoAndSwaps ;

}

Handle<DefaultProbabilityTermStructure> get_risky_curve(const Date& tradeDate,
                                                        const std::vector<Real>& cds_node_vals,
                                                        bool quoting_spreads,
                                                        Size settlement_days,
                                                        const std::vector<Period>& tenors,
                                                        const Calendar& curve_instr_calendar,
                                                        const Frequency curve_instr_freq,
                                                        const BusinessDayConvention curve_instr_conv,
                                                        const DateGeneration::Rule date_gen_rule,
                                                        const DayCounter& curve_instr_day_counter,
                                                        double recovery_rate,
                                                        const Handle<YieldTermStructure>& tsCurve,
                                                        const Real coupon
                                                        )
{
    std::vector<ext::shared_ptr<DefaultProbabilityHelper> > instruments;
    if( quoting_spreads)
    {
        for (Size i = 0; i < cds_node_vals.size(); i++) {
            instruments.push_back(ext::shared_ptr<DefaultProbabilityHelper>(new SpreadCdsHelper(
                Handle<Quote>(ext::shared_ptr<Quote>(
                    new SimpleQuote(cds_node_vals[i]))), tenors[i],
                settlement_days, curve_instr_calendar,
                curve_instr_freq, curve_instr_conv,
                date_gen_rule, curve_instr_day_counter,
                recovery_rate, tsCurve)));
        }
    }
    else{
        for (Size i = 0; i < cds_node_vals.size(); i++) {
            instruments.push_back(ext::shared_ptr<DefaultProbabilityHelper>(new UpfrontCdsHelper(
                Handle<Quote>(ext::shared_ptr<Quote>(
                    new SimpleQuote(cds_node_vals[i]))),
                coupon,
                tenors[i],
                settlement_days,
                curve_instr_calendar,
                curve_instr_freq,
                curve_instr_conv,
                date_gen_rule,
                curve_instr_day_counter,
                recovery_rate,
                tsCurve)));
        }

    }

    // Bootstrap hazard rates
    ext::shared_ptr<PiecewiseDefaultCurve<HazardRate, BackwardFlat> >
        hazardRateStructure(
            new PiecewiseDefaultCurve<HazardRate, BackwardFlat>(tradeDate, instruments,
                                                                Actual365Fixed()));

    Handle<DefaultProbabilityTermStructure> probability(hazardRateStructure);
    return probability ;

}

void example_bond()
{
    Calendar calendar = TARGET();

    Date todaysDate(1, May, 2019);
    // must be a business day
    Natural settlementDays = 2;
    Date settlementDate = calendar.advance(todaysDate, settlementDays, Days);

    // Start with today, trade is today
    Settings::instance().evaluationDate() = todaysDate;

    // Define a flat discount curve
    ext::shared_ptr<Quote> flatRate(new SimpleQuote(0.01));
    Handle<YieldTermStructure> tsCurve(
        ext::make_shared<FlatForward>(todaysDate, Handle<Quote>(flatRate), Actual365Fixed()));

    Date issueDate = Date(1, March, 2019) ;
    Date maturityDate = Date(1, March, 2024);
    DayCounter bond_day_counter = ActualActual(ActualActual::Bond) ;

    // Another discount curve using zero coupon bonds
    std::vector<Real> spot_rates = {0., 0.04, 0.04, 0.05, 0.06};
    std::vector<Date> spot_dates = {issueDate, issueDate + Period(6, Months),
                                    issueDate + Period(12, Months),
                                    issueDate + Period(2, Years),
                                    issueDate + Period(3, Years)};

    ZeroCurve dc1(spot_dates, spot_rates, bond_day_counter, calendar);

    // and one construction for the ISDA curve
    Integer isda_settlement_days = 3 ;
    Calendar dep_calendar = TARGET();
    BusinessDayConvention dep_convention = ModifiedFollowing;
    bool end_of_month = true;
    DayCounter dep_day_counter = Actual360();

    // First the deposit part:
    Rate dep1m_mkt = 0.01 ;
    Rate dep3m_mkt = 0.01 ;
    Rate dep6m_mkt = 0.01;
    Rate dep12m_mkt = 0.04;

    ext::shared_ptr<Quote> dep1m_rate(new SimpleQuote(dep1m_mkt));
    ext::shared_ptr<Quote> dep3m_rate(new SimpleQuote(dep3m_mkt));
    ext::shared_ptr<Quote> dep6m_rate(new SimpleQuote(dep6m_mkt));
    ext::shared_ptr<Quote> dep12m_rate(new SimpleQuote(dep12m_mkt));

    auto xx = 1 * Months ;
    ext::shared_ptr<RateHelper> dep1m(new DepositRateHelper(
        Handle<Quote>(dep1m_rate),
        1 * Months, isda_settlement_days,
        dep_calendar,
        dep_convention,
        end_of_month, dep_day_counter) );

    ext::shared_ptr<RateHelper> dep3m(new DepositRateHelper(
        Handle<Quote>(dep3m_rate),
        3 * Months, isda_settlement_days,
        dep_calendar,
        dep_convention,
        end_of_month, dep_day_counter) );

    ext::shared_ptr<RateHelper> dep6m(new DepositRateHelper(
        Handle<Quote>(dep6m_rate),
        6 * Months, isda_settlement_days,
        dep_calendar,
        dep_convention,
        end_of_month, dep_day_counter) );

    ext::shared_ptr<RateHelper> dep12m(new DepositRateHelper(
        Handle<Quote>(dep12m_rate),
        12 * Months, isda_settlement_days,
        dep_calendar, dep_convention,
        end_of_month, dep_day_counter) );

    // Now create the swap insturment helpers:

    Rate sw2y_mkt = 0.025;
    Rate sw3y_mkt = 0.03;
    Rate sw5y_mkt = 0.04;
    Rate sw10y_mkt = 0.043;
    Rate sw15y_mkt = 0.051;

    ext::shared_ptr<Quote> sw2y_rate(new SimpleQuote(sw2y_mkt));
    ext::shared_ptr<Quote> sw3y_rate(new SimpleQuote(sw3y_mkt));
    ext::shared_ptr<Quote> sw5y_rate(new SimpleQuote(sw5y_mkt));
    ext::shared_ptr<Quote> sw10y_rate(new SimpleQuote(sw10y_mkt));
    ext::shared_ptr<Quote> sw15y_rate(new SimpleQuote(sw15y_mkt));

    Calendar swCalendar = calendar ;
    Frequency swFixedLegFrequency = Annual;
    BusinessDayConvention swFixedLegConvention = Unadjusted;
    DayCounter swFixLegDayCounter = Thirty360(Thirty360::European);
    ext::shared_ptr<IborIndex> swFloatingLegIndex(new Euribor6M) ;

    const Period forwardStart(1*Days);

    ext::shared_ptr<RateHelper> sw2y(new SwapRateHelper(
        Handle<Quote>(sw2y_rate),
        2 * Years,
        swCalendar,
        swFixedLegFrequency,
        swFixedLegConvention,
        swFixLegDayCounter,
        swFloatingLegIndex,
        Handle<Quote>(),
        forwardStart));

    ext::shared_ptr<RateHelper> sw3y(new SwapRateHelper(
        Handle<Quote>(sw3y_rate),
        3 * Years,
        swCalendar,
        swFixedLegFrequency,
        swFixedLegConvention,
        swFixLegDayCounter,
        swFloatingLegIndex,
        Handle<Quote>(),
        forwardStart));
    ext::shared_ptr<RateHelper> sw5y(new SwapRateHelper(
        Handle<Quote>(sw5y_rate),
        5 * Years,
        swCalendar,
        swFixedLegFrequency,
        swFixedLegConvention,
        swFixLegDayCounter,
        swFloatingLegIndex,
        Handle<Quote>(),
        forwardStart));
    ext::shared_ptr<RateHelper> sw10y(new SwapRateHelper(
        Handle<Quote>(sw10y_rate),
        10 * Years,
        swCalendar,
        swFixedLegFrequency,
        swFixedLegConvention,
        swFixLegDayCounter,
        swFloatingLegIndex,
        Handle<Quote>(),
        forwardStart));

    ext::shared_ptr<RateHelper> sw15y(new SwapRateHelper(
        Handle<Quote>(sw15y_rate),
        15 * Years,
        swCalendar,
        swFixedLegFrequency,
        swFixedLegConvention,
        swFixLegDayCounter,
        swFloatingLegIndex,
        Handle<Quote>(),
        forwardStart));

    // Build the ISDA curve out of the depo and swap rates
    std::vector<ext::shared_ptr<RateHelper>> depoAndSwaps ;
    depoAndSwaps.push_back(dep1m) ;
    depoAndSwaps.push_back(dep3m) ;
    depoAndSwaps.push_back(dep6m) ;
    depoAndSwaps.push_back(dep12m) ;

    depoAndSwaps.push_back(sw2y) ;
    depoAndSwaps.push_back(sw3y) ;
    depoAndSwaps.push_back(sw5y) ;
    depoAndSwaps.push_back(sw10y) ;
    depoAndSwaps.push_back(sw15y) ;

    DayCounter termStructureDayCounter = ActualActual(ActualActual::ISDA);

    ext::shared_ptr<YieldTermStructure> isda_curve(
        new PiecewiseYieldCurve<Discount, LogLinear>(
            settlementDate, depoAndSwaps, termStructureDayCounter)
            );

    Handle<YieldTermStructure> isda_curve_handle(isda_curve) ;
    Schedule fixedBondSchedule(issueDate,
                               maturityDate,
                               Period(Quarterly),
                               calendar,
                               Unadjusted, Unadjusted,
                               DateGeneration::Backward, false);

    Real faceAmount = 100;

    FixedRateBond fixedRateBond(
        settlementDays,
        faceAmount,
        fixedBondSchedule,
        std::vector<Rate>(1, 0.05), bond_day_counter,
        ModifiedFollowing,
        100.0, Date(15, May, 2007));

    RelinkableHandle<YieldTermStructure> dcTermStructure;
    dcTermStructure.linkTo(isda_curve) ;


    ext::shared_ptr<PricingEngine> bondPricingEngine(
        new DiscountingBondEngine(dcTermStructure)) ;

    fixedRateBond.setPricingEngine(bondPricingEngine) ;
    Real price1 = fixedRateBond.NPV() ;
    Real price2 = fixedRateBond.cleanPrice() ;
    std::cout<<price2 ;

    Real mkt_price = 99. ;

    Spread zspread = BondFunctions::zSpread(fixedRateBond,
                           mkt_price,
                           //*tsCurve,
                           isda_curve,
                           bond_day_counter,
                           Compounding::Simple,
                           Semiannual,
                           settlementDate,
                           1e-6,
                           100,
                           0.01);

    // redo the calculation with a shifted by zSpread curve and see if it works
    Handle<Quote> zspread_handle( ext::shared_ptr<Quote>(new SimpleQuote(zspread))) ;


    ext::shared_ptr<YieldTermStructure> spreaded_curve(
        new ZeroSpreadedTermStructure(
        isda_curve_handle,
        zspread_handle,
        Compounding::Simple,
        Semiannual,
        bond_day_counter
        ) );

    std::cout<<zspread<<std::endl;

    Handle<YieldTermStructure> dcTermStructure_zspread(spreaded_curve);
    ext::shared_ptr<PricingEngine> bondPricingEngine2(
        new DiscountingBondEngine(dcTermStructure_zspread));

    fixedRateBond.setPricingEngine(bondPricingEngine2);
    std::cout<<fixedRateBond.cleanPrice()<<std::endl ;

}

void example_risky_bond()
{
    Calendar calendar = TARGET();

    Date evalDate(1, May, 2019);
    // must be a business day
    Natural settlementDays = 2;
    Date settlementDate = calendar.advance(evalDate, settlementDays, Days);

    Date issueDate(1, April, 2019) ;
    Date maturityDate(31, Dec, 2024) ;

    // Start with today, trade is today
    Settings::instance().evaluationDate() = evalDate;

    // All details for building out the isda riskless discounting
    Calendar dep_calendar = TARGET();
    BusinessDayConvention dep_convention = ModifiedFollowing;
    bool end_of_month = true;
    DayCounter dep_day_counter = Actual360();

    std::vector<Real> dep_rates = { 0.01, 0.015, 0.025, 0.026};

    Calendar swCalendar = calendar ;
    Frequency swFixedLegFrequency = Annual;
    BusinessDayConvention swFixedLegConvention = Unadjusted;
    DayCounter swFixLegDayCounter = Thirty360(Thirty360::European);
    ext::shared_ptr<IborIndex> swFloatingLegIndex(new Euribor6M) ;
    std::vector<Real> sw_rates = {0.025, 0.03, 0.04, 0.045, 0.055};

    std::vector<ext::shared_ptr<RateHelper>> isda_curve_nodes = get_isda_curve_nodes(
        dep_rates, sw_rates, dep_calendar, dep_convention, dep_day_counter, false,
        swCalendar, swFixedLegFrequency, swFixedLegConvention,
        swFixLegDayCounter, swFloatingLegIndex);

    DayCounter termStructureDayCounter = ActualActual(ActualActual::ISDA);
    ext::shared_ptr<YieldTermStructure> isda_curve(
        new PiecewiseYieldCurve<Discount, LogLinear>(
            settlementDate, isda_curve_nodes, termStructureDayCounter)
    );

    Handle<YieldTermStructure> isda_curve_handle(isda_curve) ;

    Size cds_settlement_days = 0 ;
    Calendar cds_curve_calendar = TARGET();
    Frequency cds_curve_freq = Quarterly ;
    BusinessDayConvention cds_curve_conv = Following ;
    DateGeneration::Rule cds_curve_date_gen_rule = QuantLib::DateGeneration::TwentiethIMM ;
    DayCounter cds_curve_day_counter = Actual365Fixed();
    Real cds_curve_rr = 0.5 ;
    Real cds_curve_coupon = 100. ;
    Size cds_curve_settlement_days = 0 ;
    std::vector<Real> cds_spreads = {0.0150, 0.02, 0.025, 0.03 };
    std::vector<Period> cds_tenors = {3 * Months, 6 * Months, 1 * Years, 2 * Years};
    Real cds_curve_coupons = 100. ;

    Handle<DefaultProbabilityTermStructure> hr_curve = get_risky_curve(evalDate,
                                                                       cds_spreads,
                                                                       true,
                                                                       cds_curve_settlement_days,
                                                                       cds_tenors,
                                                                       cds_curve_calendar,
                                                                       cds_curve_freq,
                                                                       cds_curve_conv,
                                                                       cds_curve_date_gen_rule,
                                                                       cds_curve_day_counter,
                                                                       cds_curve_rr,
                                                                       isda_curve_handle,
                                                                       cds_curve_coupons);
    Schedule fixedBondSchedule(issueDate,
                               maturityDate,
                               Period(Quarterly),
                               calendar,
                               Unadjusted, Unadjusted,
                               DateGeneration::Backward, false);

    Real coupon = 20. ;
    DayCounter fixed_bond_day_counter = Actual360();
    BusinessDayConvention fixed_bond_conv = Following ;

    std::vector<Real> notionals = {100.};

    RiskyFixedBond risky_bond = \
        RiskyFixedBond("some_name",
                       EURCurrency(),
                       cds_curve_rr,
                       hr_curve,
                       fixedBondSchedule,
                       coupon,
                       fixed_bond_day_counter,
                       fixed_bond_conv,
                       notionals,
                       isda_curve_handle) ;

    std::vector<ext::shared_ptr<CashFlow>> res = risky_bond.cashflows() ;

}

void example_flat_curve()
{
    Rate riskFreeRate = 0.05 ;
    ext::shared_ptr<Quote> rate(new SimpleQuote(riskFreeRate));
    Calendar cal = TARGET();
    Date today = cal.adjust(Date::todaysDate());
    DayCounter dayCounter = Actual365Fixed();
    Handle<YieldTermStructure> discountCurve(
        ext::shared_ptr<YieldTermStructure>(
            new FlatForward(today, Handle<Quote>(rate), dayCounter)));

    Handle<YieldTermStructure> discountCurve2(
        ext::shared_ptr<YieldTermStructure>(
            new FlatForward(0,cal, Handle<Quote>(rate), dayCounter)));

    Date some_date = Date(1, January, 2025) ;
    Date ref_date = Date(1, June, 2020) ;
    Settings::instance().evaluationDate() = ref_date;

    //std::cout<<discountCurve->discount(some_date);
    std::cout<<discountCurve2->discount(some_date);
}



int main(int argc, char *argv[]) {
    example_flat_curve();
};

