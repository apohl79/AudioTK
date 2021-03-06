
%{
#include <ATK/Tools/SumFilter.h>
%}

namespace ATK
{
  template<class DataType>
  class SumFilter: public BaseFilter
  {
  public:
    SumFilter();
    ~SumFilter();
  };
}

%template(Int16SumFilter) ATK::SumFilter<std::int16_t>;
%template(Int32SumFilter) ATK::SumFilter<std::int32_t>;
%template(Int64SumFilter) ATK::SumFilter<std::int64_t>;
%template(FloatSumFilter) ATK::SumFilter<float>;
%template(DoubleSumFilter) ATK::SumFilter<double>;
