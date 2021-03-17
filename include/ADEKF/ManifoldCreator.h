#include <boost/preprocessor.hpp>
#include <boost/preprocessor/seq.hpp>
#include <boost/preprocessor/tuple.hpp>
#include "ADEKFUtils.h"

//////////////////////////////////////////////////////////


#define ADEKF_IF_0(macro, ...)
#define ADEKF_IF_1(macro,...) macro(__VA_ARGS__)
#define ADEKF_IF_(BOOL,macro,...) ADEKF_IF_ ## BOOL (macro,__VA_ARGS__)
#define ADEKF_IF(BOOL,macro,...) ADEKF_IF_(BOOL,macro,__VA_ARGS__)


#define ADEKF_SEQ_NOT_EMPTY(seq) BOOST_PP_BOOL(BOOST_PP_DEC(BOOST_PP_SEQ_SIZE(()seq)))
#define ADEKF_IF_SEQ_NOT_EMPTY(seq,macro,...) ADEKF_IF(ADEKF_SEQ_NOT_EMPTY(seq),macro,__VA_ARGS__)
#define TRANSFORM_IF_NOT_EMPTY(macro,data,sequence) ADEKF_IF_SEQ_NOT_EMPTY(sequence,BOOST_PP_SEQ_TRANSFORM_S,1,macro,data,sequence)
#define FOR_EACH_IF_NOT_EMPTY(macro,data,sequence)  ADEKF_IF_SEQ_NOT_EMPTY(sequence,BOOST_PP_SEQ_FOR_EACH_R,1,macro,data,sequence)

#define ADEKF_APPLY_MACRO_ON_TUPLE(r, macro, tuple) macro (tuple)

//fix for empty sequence
#define BOOST_PP_SEQ_ENUM_0
#define ADEKF_TRANSFORM_COMMA(macro, entries) BOOST_PP_SEQ_ENUM(TRANSFORM_IF_NOT_EMPTY( ADEKF_APPLY_MACRO_ON_TUPLE, macro, entries))
#define ADEKF_TRANSFORM_COMMA_UNCHECKED(macro,sequence) BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM_S(1,ADEKF_APPLY_MACRO_ON_TUPLE,macro,sequence))
#define ADEKF_TRANSFORM_UNCHECKED(macro,sequence) BOOST_PP_SEQ_FOR_EACH_R(1,ADEKF_APPLY_MACRO_ON_TUPLE,macro,sequence)

#define ADEKF_TRANSFORM(macro, entries) FOR_EACH_IF_NOT_EMPTY(ADEKF_APPLY_MACRO_ON_TUPLE, macro, entries)

#define MANCREATOR_GETGLOBAL(name)  decltype(name)::GLOBAL_SIZE
#define ADEKF_RECURSIVE_(seq, output_function, step) \
BOOST_PP_FOR_1( \
		( \
				BOOST_PP_SEQ_SIZE(seq), \
				BOOST_PP_SEQ_HEAD(seq), \
				BOOST_PP_SEQ_TAIL(seq) (~), \
				0 ),\
		ADEKF_RECURSIVE_PREDICATE, step, output_function)

#define ADEKF_RECURSIVE_PREDICATE(r, state) BOOST_PP_TUPLE_ELEM(0,state)

#define ADEKF_RECURSIVE_STEP(r, state) ADEKF_RECURSIVE_STEP_IMPL state
#define ADEKF_RECURSIVE_STEP_IMPL(len, head, seq, dof) (BOOST_PP_DEC(len),BOOST_PP_SEQ_HEAD(seq),BOOST_PP_SEQ_TAIL(seq),dof + ADEKF_GETDOF(head))

#define ADEKF_RECURSIVE_STEP_GLOB(r, state) ADEKF_RECURSIVE_STEP_GLOB_IMPL state
#define ADEKF_RECURSIVE_STEP_GLOB_IMPL(len, head, seq, glob) (BOOST_PP_DEC(len),BOOST_PP_SEQ_HEAD(seq),BOOST_PP_SEQ_TAIL(seq),glob + MANCREATOR_GETGLOBAL(head))


#define ADEKF_CONSTRUCTOR_ARGS_NO_DEFAULT(name)  const typename adekf::StateInfo<decltype(name)>::type & arg_##name



#define ADEKF_RECURSIVE(seq,output_function) ADEKF_RECURSIVE_(seq,output_function,ADEKF_RECURSIVE_STEP)

#define ADEKF_RECURSIVE_COMMA_(seq,output_function,step) BOOST_PP_SEQ_ENUM(ADEKF_RECURSIVE_(seq,output_function,step))
#define ADEKF_RECURSIVE_COMMA(seq,output_function) ADEKF_RECURSIVE_COMMA_(seq,output_function,ADEKF_RECURSIVE_STEP)

/**
 * Struct to get default values differently for manifolds
 */
template<typename T>
struct defaultValue{
    inline static const T default_value{};};

/**
 * Struct to get default value (Eigen::Zero()) for Eigen Matrices
 */
template <typename DERIVED>
struct defaultValue<Eigen::MatrixBase<DERIVED> >{
    inline static const DERIVED default_value=DERIVED::Zero();
};



#define ADEKF_TEMPLATE_TYPES(name) typename TYPE_ ##name
#define ADEKF_TEMPLATE_ARGS(name) TYPE_ ##name  arg_##name

#define ADEKF_CONSTRUCTOR_SETTERS(name) name(arg_##name)
#define ADEKF_SETTERS(name) name=arg_##name;
#define ADEKF_DEFAULT_CONSTRUCTOR_SETTERS(name) name(defaultValue<typename adekf::StateInfo<decltype(name)>::type>::default_value)
#define ADEKF_COPY_CONSTRUCTOR(name) name(other.name)
#define ADEKF_ADD_DOF(name) + ADEKF_GETDOF(name)
#define ADEKF_ADD_GLOBAL_SIZE(name) + ADEKF_GETGLOBAL(name)
#define ADEKF_ASSIGN(name) name=other.name;
#define ADEKF_ARG_NAME(name)  arg_##name
#define ADEKF_STREAM_ASSIGN_VECTOR(v_members) vector_part << ADEKF_TRANSFORM_COMMA_UNCHECKED(ADEKF_ARG_NAME,v_members);

#define ADEKF_PLUS_OUTPUT(r, state) ADEKF_PLUS_OUTPUT_IMPL state
#define ADEKF_PLUS_OUTPUT_IMPL(len, head, seq, dof) (head + __delta.template segment<ADEKF_GETDOF(head)>(dof))

#define ADEKF_MINUS_OUTPUT(r, state) ADEKF_MINUS_OUTPUT_IMPL state
#define ADEKF_MINUS_OUTPUT_IMPL(len, head, seq, dof)  (head - __other.head)


#define ADEKF_OUTSTREAM_IMPL(name) << __state.name << std::endl
#define ADEKF_INSTREAM_IMPL(name) >> __state.name

#define ADEKF_MAP_OUTPUT(r,state) ADEKF_MAP_OUTPUT_IMPL state
#define ADEKF_MAP_OUTPUT_IMPL(len, head, seq, dof) , head(vector_part.template segment<ADEKF_GETDOF(head)>(dof))


#define ADEKF_CERES_PLUS(r,state) ADEKF_CERES_PLUS_IMPL state
#define ADEKF_CERES_PLUS_IMPL(len, head, seq, dof) result&=head(&x[glob],&delta[dof],&x_plus_delta[glob]); glob+=ADEKF_GETGLOBAL(head);
#define ADEKF_FROM_POINTER(r,state) ADEKF_FROM_POINTER_IMPL state
#define ADEKF_FROM_POINTER_IMPL(len, head, seq, glob) (decltype(head){&src[glob]})
#define ADEKF_TO_POINTER(r,state) ADEKF_TO_POINTER_IMPL state
#define ADEKF_TO_POINTER_IMPL(len, head, seq, glob) head.toPointer(&dest[glob]);


///////////////////////////////////////////////////////////

#define ADEKF_CONSTRUCTOR(name, m_members, v_members) \
explicit name<T>():ADEKF_TRANSFORM_COMMA(ADEKF_DEFAULT_CONSTRUCTOR_SETTERS,m_members)  BOOST_PP_COMMA_IF(ADEKF_SEQ_NOT_EMPTY(m_members)) vector_part(vector_part.Zero()) ADEKF_IF_SEQ_NOT_EMPTY(v_members,ADEKF_RECURSIVE,v_members,ADEKF_MAP_OUTPUT) {}\
name<T>( \
    ADEKF_TRANSFORM_COMMA(ADEKF_CONSTRUCTOR_ARGS_NO_DEFAULT,m_members) BOOST_PP_COMMA_IF(BOOST_PP_BITAND(ADEKF_SEQ_NOT_EMPTY(m_members),ADEKF_SEQ_NOT_EMPTY(v_members))) ADEKF_TRANSFORM_COMMA(ADEKF_CONSTRUCTOR_ARGS_NO_DEFAULT,v_members) \
    ) : \
    ADEKF_TRANSFORM_COMMA(ADEKF_CONSTRUCTOR_SETTERS,m_members)  BOOST_PP_COMMA_IF(ADEKF_SEQ_NOT_EMPTY(m_members)) vector_part() ADEKF_IF_SEQ_NOT_EMPTY(v_members,ADEKF_RECURSIVE,v_members,ADEKF_MAP_OUTPUT)  {\
     ADEKF_IF_SEQ_NOT_EMPTY(v_members,ADEKF_STREAM_ASSIGN_VECTOR,v_members)\
}\
name<T>(const name<T> & other):ADEKF_TRANSFORM_COMMA(ADEKF_COPY_CONSTRUCTOR,m_members)BOOST_PP_COMMA_IF(ADEKF_SEQ_NOT_EMPTY(m_members)) vector_part(other.vector_part) ADEKF_IF_SEQ_NOT_EMPTY(v_members,ADEKF_RECURSIVE,v_members,ADEKF_MAP_OUTPUT){}

#define ADEKF_VECTOR_CONSTRUCTOR(name, m_members, v_members)\
explicit name<T>( \
    ADEKF_TRANSFORM_COMMA(ADEKF_CONSTRUCTOR_ARGS_NO_DEFAULT,m_members) BOOST_PP_COMMA_IF(ADEKF_SEQ_NOT_EMPTY(m_members)) const decltype(vector_part) & arg_vector_part \
    ) : \
    ADEKF_TRANSFORM_COMMA(ADEKF_CONSTRUCTOR_SETTERS,m_members)  BOOST_PP_COMMA_IF(ADEKF_SEQ_NOT_EMPTY(m_members)) vector_part(arg_vector_part) ADEKF_IF_SEQ_NOT_EMPTY(v_members,ADEKF_RECURSIVE,v_members,ADEKF_MAP_OUTPUT)  {}


#define ADEKF_DEDUCTION_RESULT_IMPL(name, member) name<typename adekf::StateInfo<TYPE_##member>::ScalarType >
#define ADEKF_DEDUCTION_RESULT(name, member) ADEKF_DEDUCTION_RESULT_IMPL(name,member)

#define ADEKF_DEDUCTION_GUIDE(name, members) \
        template< ADEKF_TRANSFORM_COMMA(ADEKF_TEMPLATE_TYPES,members) >\
        name(ADEKF_TRANSFORM_COMMA(ADEKF_TEMPLATE_ARGS,members)) -> ADEKF_DEDUCTION_RESULT(name,BOOST_PP_SEQ_HEAD(members));


#define ADEKF_ASSIGN_OP(name, m_members)\
void operator=(const name<T> & other){\
  ADEKF_TRANSFORM(ADEKF_ASSIGN,m_members)\
  vector_part=other.vector_part;\
}\
void operator=(const name<T> && other){\
  ADEKF_TRANSFORM(ADEKF_ASSIGN,m_members)\
  vector_part=other.vector_part;\
}


#define VECTOR_PART_BOXPLUS(m_members) BOOST_PP_COMMA_IF(ADEKF_SEQ_NOT_EMPTY(m_members)) vector_part+__delta.template tail<VEC_DOF>()

#define ADEKF_BOXPLUS(name, m_members, v_members) \
template<typename Derived> \
name<ADEKF_PLUSRESULT(T,typename Derived::Scalar)> operator+(const Eigen::MatrixBase<Derived>& __delta) const { \
EIGEN_STATIC_ASSERT_SAME_VECTOR_SIZE(Derived, DeltaType<T>)\
     return name<ADEKF_PLUSRESULT(T,typename Derived::Scalar)>{ADEKF_IF_SEQ_NOT_EMPTY(m_members,ADEKF_RECURSIVE_COMMA,m_members,ADEKF_PLUS_OUTPUT) ADEKF_IF_SEQ_NOT_EMPTY(v_members,VECTOR_PART_BOXPLUS,m_members) }; \
}


#define ADEKF_BOXMINUS(name, m_members, v_members) \
template<typename T2> \
DeltaType<ADEKF_MINUSRESULT(T,T2)> operator-(const name<T2>& __other) const { \
    DeltaType<ADEKF_MINUSRESULT(T,T2)> __result; __result << ADEKF_IF_SEQ_NOT_EMPTY(m_members,ADEKF_RECURSIVE_COMMA,m_members,ADEKF_MINUS_OUTPUT) BOOST_PP_COMMA_IF(ADEKF_SEQ_NOT_EMPTY(m_members)) vector_part-__other.vector_part ;return __result;\
}

#define CALL_FUNCTION_ON_MEMBER(r, data, name)  BOOST_PP_TUPLE_ELEM(0,data) (name,BOOST_PP_TUPLE_ELEM(1,data)...);
#define CALL_FUNCTION_ON_MEMBER_WITH_OTHER(r, data, name)  BOOST_PP_TUPLE_ELEM(0,data) (name,BOOST_PP_TUPLE_ELEM(2,data). name,BOOST_PP_TUPLE_ELEM(1,data)...);

#define ADEKF_FOREACH(name, m_members, function_name)  \
template<typename FUNCTOR, typename ... ARGS>\
void function_name(FUNCTOR & functor,  ARGS &&... args) const{\
    FOR_EACH_IF_NOT_EMPTY(CALL_FUNCTION_ON_MEMBER,(functor,args),m_members) \
}\
template<typename FUNCTOR, typename OtherScalar, typename ... ARGS>\
void function_name ##WithOther (FUNCTOR & functor,const name<OtherScalar> & other ,ARGS && ... args) const {\
FOR_EACH_IF_NOT_EMPTY(CALL_FUNCTION_ON_MEMBER_WITH_OTHER,(functor,args,other),m_members)\
}


#define ADEKF_OUTSTREAM(name, m_members) \
friend std::ostream& operator<<(std::ostream& __stream, const name& __state) { \
    return __stream ADEKF_TRANSFORM(ADEKF_OUTSTREAM_IMPL,m_members) << __state.vector_part; \
}

#define ADEKF_INSTREAM(name, members) \
friend std::istream& operator>>(std::istream& __stream, name& __state) { \
    return __stream ADEKF_TRANSFORM(ADEKF_INSTREAM_IMPL,members); \
}
//Compatibility with ceres optimizers
#ifdef MANIFOLD_WITH_CERES
#include <ceres/local_parameterization.h>
#include <ceres/autodiff_local_parameterization.h>
#define CERES_ADAPTER(name, m_members)\
static constexpr unsigned MAN_GLOBAL_SIZE=0 ADEKF_TRANSFORM_UNCHECKED(ADEKF_ADD_GLOBAL_SIZE,m_members);\
static constexpr unsigned GLOBAL_SIZE=MAN_GLOBAL_SIZE+VEC_DOF;\
template<typename TYPE>\
bool operator()(const TYPE* x,const TYPE* delta,TYPE* x_plus_delta) const {\
    size_t glob=0;\
    bool result=true;\
    ADEKF_RECURSIVE(m_members,ADEKF_CERES_PLUS)\
    return result;\
}\
static ceres::LocalParameterization * getLocalParameterization(){ \
    return new ceres::ProductParameterization{new ceres::AutoDiffLocalParameterization<name,MAN_GLOBAL_SIZE,MAN_DOF>{}, new ceres::IdentityParameterization{VEC_DOF}};\
}\
explicit name<T>(const T* const src): name{ADEKF_RECURSIVE_COMMA_(m_members,ADEKF_FROM_POINTER,ADEKF_RECURSIVE_STEP_GLOB) BOOST_PP_COMMA_IF(ADEKF_SEQ_NOT_EMPTY(m_members))  Eigen::Map<const Eigen::Matrix<T,VEC_DOF,1>>(&src[MAN_GLOBAL_SIZE])}{}\
void toPointer(T* dest){\
  ADEKF_RECURSIVE_(m_members,ADEKF_TO_POINTER,ADEKF_RECURSIVE_STEP_GLOB)  \
  (Eigen::Map<Eigen::Matrix<T,VEC_DOF,1>>(&dest[MAN_GLOBAL_SIZE]))=vector_part;\
}
#else
#define CERES_ADAPTER(name, m_members)
#endif


////////////////////////////////////////////////////////////////////////

#define ADEKF_CONSTRUCT_MANIFOLD_INTERNAL(name, m_members, v_members) \
ADEKF_CONSTRUCTOR(name, m_members, v_members) \
BOOST_PP_REMOVE_PARENS(BOOST_PP_EXPR_IF(BOOST_PP_GREATER(BOOST_PP_SEQ_SIZE(v_members),1),(ADEKF_VECTOR_CONSTRUCTOR(name,m_members,v_members))))\
using ScalarType = T; \
static constexpr unsigned MAN_DOF=0 ADEKF_TRANSFORM(ADEKF_ADD_DOF,m_members);\
static constexpr unsigned DOF = VEC_DOF +MAN_DOF; \
template<typename DT>\
using DeltaType=Eigen::Matrix<DT,DOF,1>;\
ADEKF_BOXPLUS(name, m_members, v_members) \
ADEKF_BOXMINUS(name, m_members, v_members) \
ADEKF_ASSIGN_OP(name,m_members)\
ADEKF_FOREACH(name,m_members v_members, forEach)\
ADEKF_FOREACH(name,m_members, forEachManifold)\
ADEKF_FOREACH(name,v_members, forEachVector)\
ADEKF_OUTSTREAM(name, m_members) \
ADEKF_IF_SEQ_NOT_EMPTY(m_members,CERES_ADAPTER,name,m_members) \
EIGEN_MAKE_ALIGNED_OPERATOR_NEW\
};\
	ADEKF_DEDUCTION_GUIDE(name,m_members v_members)\
	typedef adekf::CovarianceOf<name<double>> name ##Cov;

//ADEKF_INSTREAM(name, m_members, v_members) \



#define ADEKF_APPLY_MACRO_ON_DUPLET(r, macro, tuple) macro   tuple

/**
 * Those structs are required to break down the Eigen template arguments so that you do not have to write commas in the macros
 */
template<int VEC_DOF, int SIZE>
struct MatrixBlock{
template<typename T>
using type=Eigen::Block<Eigen::Matrix<T,VEC_DOF,1>, SIZE, 1, false>;
};
	template<int  VEC_DOF>
struct MatrixBlockWrapper{
    template<int SIZE>
    using type=MatrixBlock<VEC_DOF,SIZE>;
};

#define CREATE_ATTRIBUTE(TYPE,NAME) TYPE<T> NAME;

#define CREATE_VECTOR_ATTRIBUTE(SIZE,NAME) MatrixBlockWrapper<VEC_DOF>::template type<SIZE>::template type<T> NAME;

#define RETURN_NAME(TYPE, NAME) (NAME)

#define ADD_UP_SIZE(SIZE,NAME) +SIZE

#define ADEKF_UNZIP_ATTRIBUTES(macro,seq) FOR_EACH_IF_NOT_EMPTY(ADEKF_APPLY_MACRO_ON_DUPLET,macro,seq)


#define ADEKF_MANIFOLD_BOTH(name, manifolds, vectors) template<typename T> struct name :public adekf::CompoundManifold\
{\
    static constexpr unsigned VEC_DOF=0 ADEKF_UNZIP_ATTRIBUTES(ADD_UP_SIZE,vectors) ;\
    /**   \
    * Contains all vector data of this CompoundManifold  \
    */\
    ADEKF_UNZIP_ATTRIBUTES(CREATE_ATTRIBUTE,manifolds)             \
    Eigen::Matrix<T,VEC_DOF,1> vector_part; \
    ADEKF_UNZIP_ATTRIBUTES(CREATE_VECTOR_ATTRIBUTE,vectors)    \
	ADEKF_CONSTRUCT_MANIFOLD_INTERNAL(name,   ADEKF_UNZIP_ATTRIBUTES(RETURN_NAME,manifolds) , ADEKF_UNZIP_ATTRIBUTES(RETURN_NAME,vectors)) 																	\



#define EXPAND_1(name,manifolds,...) ADEKF_MANIFOLD_BOTH(name,manifolds,BOOST_PP_VARIADIC_TO_SEQ( __VA_ARGS__))
#define EXPAND_0(name,manifolds,...) ADEKF_MANIFOLD_BOTH(name,manifolds,)
#define EXPAND_(number, ... )  EXPAND_##number(__VA_ARGS__)

#define ADEKF_PARSE_MANIFOLD(number,  ...) EXPAND_(number,__VA_ARGS__)

#define ADEKF_MANIFOLD(...) ADEKF_PARSE_MANIFOLD(BOOST_PP_BOOL(BOOST_PP_DEC(BOOST_PP_DEC(BOOST_PP_VARIADIC_SIZE(__VA_ARGS__)))),__VA_ARGS__)
