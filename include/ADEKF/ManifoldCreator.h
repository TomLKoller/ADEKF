#include <boost/preprocessor.hpp>
#include <boost/preprocessor/seq.hpp>
#include <boost/preprocessor/tuple.hpp>
#include "ADEKFUtils.h"
#define ADEKF_PLUSRESULT(T1, T2) decltype(std::declval<T1>() + std::declval<T2>())

#define ADEKF_MINUSRESULT(T1, T2) decltype(std::declval<T1>() - std::declval<T2>())

//////////////////////////////////////////////////////////


#define ADEKF_IF_0(macro, ...)
#define ADEKF_IF_1(macro,...) macro(__VA_ARGS__)
#define ADEKF_IF_(BOOL,macro,...) ADEKF_IF_ ## BOOL (macro,__VA_ARGS__)
#define ADEKF_IF(BOOL,macro,...) ADEKF_IF_(BOOL,macro,__VA_ARGS__)


#define ADEKF_SEQ_NOT_EMPTY(seq) BOOST_PP_BOOL(BOOST_PP_DEC(BOOST_PP_SEQ_SIZE(()seq)))
#define ADEKF_IF_SEQ_NOT_EMPTY(seq,macro,...) ADEKF_IF(ADEKF_SEQ_NOT_EMPTY(seq),macro,__VA_ARGS__)
#define TRANSFORM_IF_NOT_EMPTY(macro,data,sequence) ADEKF_IF_SEQ_NOT_EMPTY(sequence,BOOST_PP_SEQ_TRANSFORM_S,1,macro,data,sequence)
#define FOR_EACH_IF_NOT_EMPTY(macro,data,sequence)  ADEKF_IF_SEQ_NOT_EMPTY(sequence,BOOST_PP_SEQ_FOR_EACH_R,1,macro,data,sequence)

#define ADEKF_TO_ENUM(r,macro,data,sequence) BOOST_PP_SEQ_ENUM(sequence)

#define ADEKF_APPLY_MACRO_ON_TUPLE(r, macro, tuple) macro (tuple)

//fix for empty sequence
#define BOOST_PP_SEQ_ENUM_0
#define ADEKF_TRANSFORM_COMMA(macro, entries) BOOST_PP_SEQ_ENUM(TRANSFORM_IF_NOT_EMPTY( ADEKF_APPLY_MACRO_ON_TUPLE, macro, entries))

#define ADEKF_TRANSFORM(macro, entries) FOR_EACH_IF_NOT_EMPTY(ADEKF_APPLY_MACRO_ON_TUPLE, macro, entries)



#define ADEKF_RECURSIVE(seq, output_function) \
BOOST_PP_FOR_1( \
		( \
				BOOST_PP_SEQ_SIZE(seq), \
				BOOST_PP_SEQ_HEAD(seq), \
				BOOST_PP_SEQ_TAIL(seq) (~), \
				0 ),\
		ADEKF_RECURSIVE_PREDICATE, ADEKF_RECURSIVE_STEP, output_function)

#define ADEKF_RECURSIVE_PREDICATE(r, state) BOOST_PP_TUPLE_ELEM(0,state)

#define ADEKF_RECURSIVE_STEP(r, state) ADEKF_RECURSIVE_STEP_IMPL state
#define ADEKF_RECURSIVE_STEP_IMPL(len, head, seq, dof) (BOOST_PP_DEC(len),BOOST_PP_SEQ_HEAD(seq),BOOST_PP_SEQ_TAIL(seq),dof + ADEKF_GETDOF(head))

#define ADEKF_CONSTRUCTOR_ARGS_NO_DEFAULT(name)  const typename adekf::StateInfo<decltype(name)>::type & arg_##name


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


#define ADEKF_CONSTRUCTOR_ARGS(name) ADEKF_CONSTRUCTOR_ARGS_NO_DEFAULT(name)=defaultValue<typename adekf::StateInfo<decltype(name)>::type>::default_value

#define ADEKF_TEMPLATE_TYPES(name) typename TYPE_ ##name
#define ADEKF_TEMPLATE_ARGS(name) TYPE_ ##name  arg_##name

#define ADEKF_CONSTRUCTOR_SETTERS(name) name(arg_##name)
#define ADEKF_SETTERS(name) name=arg_##name;
#define ADEKF_DEFAULT_SETTERS(name) name=defaultValue<typename adekf::StateInfo<decltype(name)>::type>::default_value;
#define ADEKF_DEFAULT_CONSTRUCTOR_SETTERS(name) name(defaultValue<typename adekf::StateInfo<decltype(name)>::type>::default_value)
#define ADEKF_COPY_CONSTRUCTOR(name) name(other.name)
#define ADEKF_ADD_DOF(name) + ADEKF_GETDOF(name)

#define ADEKF_PLUS_OUTPUT(r, state) ADEKF_PLUS_OUTPUT_IMPL state
#define ADEKF_PLUS_OUTPUT_IMPL(len, head, seq, dof) (head + __delta.template segment<ADEKF_GETDOF(head)>(dof))

#define ADEKF_MINUS_OUTPUT(r, state) ADEKF_MINUS_OUTPUT_IMPL state
#define ADEKF_MINUS_OUTPUT_IMPL(len, head, seq, dof) __result.template segment<ADEKF_GETDOF(head)>(dof) = head - __other.head;


#define ADEKF_OUTSTREAM_IMPL(name) << __state.name << std::endl
#define ADEKF_INSTREAM_IMPL(name) >> __state.name

#define ADEKF_MAP_OUTPUT(r,state) ADEKF_MAP_OUTPUT_IMPL state
#define ADEKF_MAP_OUTPUT_IMPL(len, head, seq, dof) , head(vector_part.template block<ADEKF_GETDOF(head),1>(dof,0))

///////////////////////////////////////////////////////////

#define ADEKF_CONSTRUCTOR(name, m_members,v_members) \
name<T>():ADEKF_TRANSFORM_COMMA(ADEKF_DEFAULT_CONSTRUCTOR_SETTERS,m_members)  BOOST_PP_COMMA_IF(ADEKF_SEQ_NOT_EMPTY(m_members)) vector_part() ADEKF_IF_SEQ_NOT_EMPTY(v_members,ADEKF_RECURSIVE,v_members,ADEKF_MAP_OUTPUT) {\
ADEKF_TRANSFORM(ADEKF_DEFAULT_SETTERS,v_members)}\
name<T>( \
    ADEKF_TRANSFORM_COMMA(ADEKF_CONSTRUCTOR_ARGS_NO_DEFAULT,m_members) BOOST_PP_COMMA_IF(BOOST_PP_BITAND(ADEKF_SEQ_NOT_EMPTY(m_members),ADEKF_SEQ_NOT_EMPTY(v_members))) ADEKF_TRANSFORM_COMMA(ADEKF_CONSTRUCTOR_ARGS_NO_DEFAULT,v_members) \
    ) : \
    ADEKF_TRANSFORM_COMMA(ADEKF_CONSTRUCTOR_SETTERS,m_members)  BOOST_PP_COMMA_IF(ADEKF_SEQ_NOT_EMPTY(m_members)) vector_part() ADEKF_IF_SEQ_NOT_EMPTY(v_members,ADEKF_RECURSIVE,v_members,ADEKF_MAP_OUTPUT)  {\
    ADEKF_TRANSFORM(ADEKF_SETTERS,v_members)\
}\
name<T>(const name<T> & other):ADEKF_TRANSFORM_COMMA(ADEKF_COPY_CONSTRUCTOR,m_members)BOOST_PP_COMMA_IF(ADEKF_SEQ_NOT_EMPTY(m_members)) vector_part(other.vector_part) ADEKF_IF_SEQ_NOT_EMPTY(v_members,ADEKF_RECURSIVE,v_members,ADEKF_MAP_OUTPUT){}

#define ADEKF_VECTOR_CONSTRUCTOR(name, m_members,v_members)\
name<T>( \
    ADEKF_TRANSFORM_COMMA(ADEKF_CONSTRUCTOR_ARGS_NO_DEFAULT,m_members) BOOST_PP_COMMA_IF(ADEKF_SEQ_NOT_EMPTY(m_members)) const decltype(vector_part) & arg_vector_part \
    ) : \
    ADEKF_TRANSFORM_COMMA(ADEKF_CONSTRUCTOR_SETTERS,m_members)  BOOST_PP_COMMA_IF(ADEKF_SEQ_NOT_EMPTY(m_members)) vector_part(arg_vector_part) ADEKF_IF_SEQ_NOT_EMPTY(v_members,ADEKF_RECURSIVE,v_members,ADEKF_MAP_OUTPUT)  {}


#define ADEKF_DEDUCTION_RESULT_IMPL(name,member) name<typename adekf::StateInfo<TYPE_##member>::ScalarType >
#define ADEKF_DEDUCTION_RESULT(name,member) ADEKF_DEDUCTION_RESULT_IMPL(name,member)

#define ADEKF_DEDUCTION_GUIDE(name, members) \
		template< ADEKF_TRANSFORM_COMMA(ADEKF_TEMPLATE_TYPES,members) >\
		name(ADEKF_TRANSFORM_COMMA(ADEKF_TEMPLATE_ARGS,members)) -> ADEKF_DEDUCTION_RESULT(name,BOOST_PP_SEQ_HEAD(members));




/*#define ADEKF_BOXPLUS(name, members) \
template<typename T2> \
name<ADEKF_PLUSRESULT(T,T2)> operator+(const Eigen::Matrix<T2, DOF, 1>& __delta) const { \
    return name<ADEKF_PLUSRESULT(T,T2)>(BOOST_PP_SEQ_ENUM(ADEKF_RECURSIVE(members,ADEKF_PLUS_OUTPUT))); \
}


#define ADEKF_BOXMINUS(name, members) \
template<typename T2> \
Eigen::Matrix<ADEKF_MINUSRESULT(T,T2),DOF,1> operator-(const name<T2>& __other) const { \
    Eigen::Matrix<ADEKF_MINUSRESULT(T,T2),DOF,1> __result;ADEKF_RECURSIVE(members,ADEKF_MINUS_OUTPUT)return __result;\
}*/


#define ADEKF_BOXPLUS(name, m_members, v_members) \
template<typename T2> \
name<ADEKF_PLUSRESULT(T,T2)> operator+(const Eigen::Matrix<T2, DOF, 1>& __delta) const { \
     return name<ADEKF_PLUSRESULT(T,T2)>{ADEKF_IF_SEQ_NOT_EMPTY(m_members,ADEKF_RECURSIVE,m_members,ADEKF_PLUS_OUTPUT),vector_part+__delta.template segment<VEC_DOF>(MAN_DOF)}; \
}


#define ADEKF_BOXMINUS(name, m_members,v_members) \
template<typename T2> \
Eigen::Matrix<ADEKF_MINUSRESULT(T,T2),DOF,1> operator-(const name<T2>& __other) const { \
    Eigen::Matrix<ADEKF_MINUSRESULT(T,T2),DOF,1> __result;ADEKF_IF_SEQ_NOT_EMPTY(m_members,ADEKF_RECURSIVE,m_members,ADEKF_MINUS_OUTPUT) __result.template segment<VEC_DOF>(MAN_DOF)=vector_part-__other.vector_part;return __result;\
}

#define CALL_FUNCTION_ON_MEMBER(r,data,name)  BOOST_PP_TUPLE_ELEM(0,data) (name,BOOST_PP_TUPLE_ELEM(1,data)...);
#define CALL_FUNCTION_ON_MEMBER_WITH_OTHER(r,data,name)  BOOST_PP_TUPLE_ELEM(0,data) (name,BOOST_PP_TUPLE_ELEM(2,data). name,BOOST_PP_TUPLE_ELEM(1,data)...);

#define ADEKF_FOREACH(name,m_members, function_name)  \
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

////////////////////////////////////////////////////////////////////////

#define ADEKF_CONSTRUCT_MANIFOLD_INTERNAL(name, m_members, v_members) \
ADEKF_CONSTRUCTOR(name, m_members, v_members) \
BOOST_PP_REMOVE_PARENS(BOOST_PP_EXPR_IF(BOOST_PP_GREATER(BOOST_PP_SEQ_SIZE(v_members),1),(ADEKF_VECTOR_CONSTRUCTOR(name,m_members,v_members))))\
using ScalarType = T; \
static constexpr unsigned MAN_DOF=0 ADEKF_TRANSFORM(ADEKF_ADD_DOF,m_members);\
static constexpr unsigned DOF = VEC_DOF +MAN_DOF; \
ADEKF_BOXPLUS(name, m_members, v_members) \
ADEKF_BOXMINUS(name, m_members, v_members) \
ADEKF_FOREACH(name,m_members v_members, forEach)\
ADEKF_FOREACH(name,m_members, forEachManifold)\
ADEKF_FOREACH(name,v_members, forEachVector)\
ADEKF_OUTSTREAM(name, m_members) \
EIGEN_MAKE_ALIGNED_OPERATOR_NEW\
};\
	ADEKF_DEDUCTION_GUIDE(name,m_members v_members)

//ADEKF_INSTREAM(name, m_members, v_members) \



#define ADEKF_APPLY_MACRO_ON_DUPLET(r, macro, tuple) macro   tuple  //(BOOST_PP_TUPLE_ELEM(0,tuple),BOOST_PP_TUPLE_ELEM(1,tuple))

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
    Eigen::Matrix<T,VEC_DOF,1> vector_part; \
    ADEKF_UNZIP_ATTRIBUTES(CREATE_ATTRIBUTE,manifolds)             \
    ADEKF_UNZIP_ATTRIBUTES(CREATE_VECTOR_ATTRIBUTE,vectors)    \
	ADEKF_CONSTRUCT_MANIFOLD_INTERNAL(name,   ADEKF_UNZIP_ATTRIBUTES(RETURN_NAME,manifolds) , ADEKF_UNZIP_ATTRIBUTES(RETURN_NAME,vectors)) 																	\


#define ADEKF_STATE(name, ...) template <typename T> struct name :public adekf::CompoundManifold{\
		ADEKF_UNZIP_ATTRIBUTES(CREATE_VECTOR_ATTRIBUTE,BOOST_PP_VARIADIC_TO_SEQ( __VA_ARGS__))	\
		ADEKF_CONSTRUCT_MANIFOLD_INTERNAL(name, ADEKF_UNZIP_ATTRIBUTES(RETURN_NAME,BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)),)


#define EXPAND_1(name,manifolds,...) ADEKF_MANIFOLD_BOTH(name,manifolds,BOOST_PP_VARIADIC_TO_SEQ( __VA_ARGS__))
#define EXPAND_0(name,manifolds,...) ADEKF_MANIFOLD_BOTH(name,manifolds,)
#define EXPAND_(number, ... )  EXPAND_##number(__VA_ARGS__)

#define ADEKF_PARSE_MANIFOLD(number,  ...) EXPAND_(number,__VA_ARGS__)

#define ADEKF_MANIFOLD(...) ADEKF_PARSE_MANIFOLD(BOOST_PP_BOOL(BOOST_PP_DEC(BOOST_PP_DEC(BOOST_PP_VARIADIC_SIZE(__VA_ARGS__)))),__VA_ARGS__)
