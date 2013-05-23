package Heap::Simple;
use Carp;
use strict;

# Switch selecting XS or pure perl
use vars qw($VERSION @ISA @implementors);
$VERSION = "0.13";

unless (@ISA) {
    @implementors = qw(Heap::Simple::XS(0.10) Heap::Simple::Perl(0.13))
        unless @implementors;
    for my $i (@implementors) {
        my $plugin = $i;
        my $version = $plugin =~ s/\((.+)\)\z// ? $1 : undef;
        my $p = "$plugin.pm";
        $p =~ s!::!/!g;
        if (eval { require $p; 1 }) {
            *new = $plugin->can("new") ||
                croak "Package '$plugin' does not support 'new'";
            if (defined $version) {
                if (defined(my $v = $plugin->VERSION)) {
                    if ($v < $version) {
                        carp "Need '$plugin' version $version, found version $v";
                        next;
                    }
                } else {
                    croak "Could not determine the version of '$plugin'";
                }
            }
            @ISA = ($plugin);
            last;
        }
    }
    @ISA || croak "Can't locate ", join(" nor ", @implementors), " in \@INC. This probably means you didn't install any of them. (\@INC contains: @INC)";
}

1;
__END__

=head1 NAME

Heap::Simple - Fast and easy to use classic heaps

=head1 SYNOPSIS

    use Heap::Simple;

    # Create a heap
    my $heap = Heap::Simple->new;
    my $heap = Heap::Simple->new(%options);

    # Put data in the heap
    $heap->insert(@new_elements);
    # Put data in a "Object" or "Any" heap with a given key
    $heap->key_insert($key1, $element1, $key2, $element2, ...);

    # Extract the top value
    $element = $heap->extract_top;	# croaks on an empty heap
    $element = $heap->extract_first;	# returns undef on an empty heap

    # Get the top value but leave it in the heap
    $element = $heap->top;		# croaks on an empty heap
    $element = $heap->first;		# returns undef on an empty heap

    # Find the top key in the heap
    $top_key = $heap->top_key;	  	# return infinity on an empty heap
					# croaks if there's no infinity
    $top_key = $heap->first_key;  	# returns undef   on an empty heap

    # Ordered extract of all data whose key is not above a given value
    @elements = $heap->extract_upto($max_key);

    # Ordered extract of all data
    @elements = $heap->extract_all;

    # Empty the heap
    $heap->clear;

    # Find the number of elements
    $count = $heap->count;

    # Get all keys (not sorted)
    @keys = $heap->keys;
    # Get all values (not sorted)
    @values = $heap->values;

    # Find the key corresponding to a value
    $key = $heap->key($value);

    # Get/Set user_data
    $user_data  = $heap->user_data;
    $old_data   = $heap->user_data($new_data);

    # Get/Set infinity
    $infinity     = $heap->infinity;
    $old_infinity = $heap->infinity($new_data);

    # Get the position of a key in an element
    $key_index    = $heap->key_index;
    $key_name     = $heap->key_name;
    $key_method   = $heap->key_method;
    $key_function = $heap->key_function;

    # Return the value of several things that were set in new:
    $wrapped      = $heap->wrapped;
    $max_count    = $heap->max_count;
    $can_die      = $heap->can_die;
    $dirty        = $heap->dirty;
    $order	  = $heap->order;
    @elements     = $heap->elements;
    $elements     = $heap->elements;

    # Move all elements out of each heap in @heaps and into $heap
    $heap->absorb(@heaps);	# As if doing a repeated $heap->insert
    $heap->key_absorb(@heaps);	# As if doing a repeated $heap->key_insert

    # merge already sorted arrays into a new sorted array
    # This doesn't disturb the elements already in the heap
    my $merged_aref = $heap->merge_arrays($aref1, $aref2, ...);

    # Which class does the actual work ?
    $implementation = Heap::Simple->implementation;

=head1 EXAMPLE1

When key and value are kept separate:

    use Heap::Simple;
    my $heap = Heap::Simple->new(elements => "Any");

    $heap->key_insert(8, "bar");
    $heap->key_insert(5, "foo");

    # This will print foo (5 is the lowest key)
    print "First value is ", $heap->extract_top, "\n";

    $heap->key_insert(7, "baz");

    # This will print baz (7 is the lowest key)
    print "Next value is ", $heap->extract_top, "\n";
    # This will print bar (8 is now the lowest key)
    print "Next value is ", $heap->extract_top, "\n";

=head1 EXAMPLE2

When the key is part of the value:

    # This is purely for display, ignore it
    use Data::Dumper;
    $Data::Dumper::Indent = 0;
    $Data::Dumper::Terse = 1;

    # Real code starts here
    use Heap::Simple;
    my $heap = Heap::Simple->new(elements => "Array");

    $heap->insert([8, "bar"]);
    $heap->insert([5, "foo"]);

    # This will print [5, foo] (5 is the lowest key)
    print "First value is ", Dumper($heap->extract_top), "\n";

    $heap->insert([7, "baz"]);

    # This will print [7, baz] (7 is the lowest key)
    print "Next value is ", Dumper($heap->extract_top), "\n";
    # This will print [8, bar] (8 is now the lowest key)
    print "Next value is ", Dumper($heap->extract_top), "\n";

=head1 DESCRIPTION

A heap is a partially sorted structure where it's always easy to extract the
smallest element. If the collection of elements is changing dynamically, a
heap has less overhead than keeping the collection fully sorted.

The order in which equal elements get extracted is unspecified.

The main order relations supported by this module are "<" (numeric compare)
and "lt" (string compare).

The module allows you to manage data where the elements are of several
allowed types, in particular array references, hash references, objects
or just the keys themselves.

So L<new|"new"> has a lot of ways to specify element types, but the right
choices follows quite directly from the data you'll put in the heap. If the key
is part of the data (or easily derived from the data), choose an element
type that tells how to get the key out of the data, and insert elements
using L<insert|/"insert">. If the key is independent from the data or
you want to avoid repeated key calculations, use the L<Any|/"Any"> element
type and insert elements using L<key_insert|/"key_insert">.

The internals of the module do nothing with the elements inserted except
inspecting the key. This means that if you for example store a blessed
object, that's what you will get back on extract. It's also ok to keep
references to the elements around and make changes to them while they are
in the heap as long as you don't change the key.

Heap::Simple itself is just a loader for the code that will actually implement
the functionality mentioned above. You will need to install something like
L<Heap::Simple::XS|Heap::Simple::XS> or
L<Heap::Simple::Perl|Heap::Simple::Perl> to be able to actually do anything.

=head1 EXPORT

None.

=head1 METHODS

All methods that can fail will thrown an exception in case of failure
unless otherwise specified. For example, you don't have to explicitely
check the result of L<new|"new">, it will already thrown an exception in
case of bad arguments.

=over

=item X<new>my $heap = Heap::Simple->new

This simplest form creates a new (empty) heap object able to hold numeric keys.

You could for example use this to print a list of numbers from low to high:

    use Heap::Simple;

    my $heap = Heap::Simple->new;
    $heap->insert(8, 3, 14, -1, 3);
    print $heap->extract_top, " " for 1..$heap->count;
    print "\n";
    # Will print: -1 3 3 8 14

This example is silly of course. You could just as well directly use
L<perl sort|perlfunc/"sort">. But in real applications you would do the
inserting interleaved with extracting and always keeping the list sorted
would become inefficient for big lists. That is where you would use a heap.
The examples we give will however be like the one above so you can quickly
see the way in which the methods are supposed to be called.

For some applications this basic usage where you just store numeric keys will
be good enough, but usually you want to be able to store more complex elements.

Several options can help you with that:

=over 2

=item X<new_order>order => $order

$order indicates what operation is used to compare keys. Supported orders are:

=over 2

=item E<lt>

Indicates that keys are compared as numbers, and extraction is lowest value
first. This is actually the default order, so the example above could have
used:

    my $heap = Heap::Simple->new(order => "<");

and the result would have been exactly the same.

The default infinity for this order is +inf.

=item E<gt>

Indicates that keys are compared as numbers, and extraction is highest value
first.

Repeating the example with this order gives:

    use Heap::Simple;

    my $heap = Heap::Simple->new(order => ">");
    $heap->insert(8, 3, 14, -1, 3);
    print $heap->extract_top, " " for 1..$heap->count;
    print "\n";
    # Will print: 14 8 3 3 -1

The default infinity for this order is -inf.

=item lt

Indicates that the keys are compared as strings, and extraction is lowest
value first. So we could modify the "<" example to:

    use Heap::Simple;

    my $heap = Heap::Simple->new(order => "lt");
    $heap->insert("ate", 8, 3, "zzzz", 14, -1, 3, "at");
    print $heap->extract_top, " " for 1..$heap->count;
    print "\n";
    # Will print: -1 14 3 3 8 at ate zzzz

Notice how 14 comes before 3 as you would expect in lexical sorting.

The default infinity for this order is "undef" (there is no maximum string)

=item gt

Indicates that the keys are compared as strings, and extraction is highest
value first. The concept of "minimum" again becomes rather confusing.
The standard example now becomes:

    use Heap::Simple;

    my $heap = Heap::Simple->new(order => "gt");
    $heap->insert("ate", 8, 3, "zzzz", 14, -1, 3, "at");
    print $heap->extract_top, " " for 1..$heap->count;
    print "\n";
    # Will print: zzzz ate at 8 3 3 14 -1

The default infinity for this order is "" (the empty string)

=item X<CODE>$code_reference

If your keys are completely weird things, ordered neither as numbers nor as
strings and you need a special compare function, you can use this general
ordering type.

Every time two keys need to be compared, the given code reference will be
called like:

    $less = $code_reference->($key1, $key2);

This should return a true value if $key1 is smaller than $key2 and a false
value otherwise. $code_reference should imply a total order relation, so it
needs to be transitive.

Since in this case nothing can be determined about the key type, there will
be no infinity by default (even if the keys are numbers).

Example:

    use Heap::Simple;

    sub more { return $_[0] > $_[1] }

    my $heap = Heap::Simple->new(order => \&more);
    $heap->insert(8, 3, 14, -1, 3);
    print $heap->extract_top, " " for 1..$heap->count;
    print "\n";
    # Will print: 14 8 3 3 -1

The code reference will be called many times during normal heap operations
(O(log n) times for a single insert or extract on a size n heap), so only use
this order type within reason. Usually it's better to precalculate some number
or string representation of some sort of key and use normal compares on these.
You can use the L<Any element type|"Any"> and L<key_insert|"key_insert"> to
wrap the precalculated key with the corresponding element, or you can delegate
the key calculation to the L<insert|"insert"> method and use one of the
L<Method|"Method">, L<Object|"Object"> or L<Function|"Function"> element types.

Here's an example of such "fake" keys:

    # "human" sorting mixed strings
    use Heap::Simple;

    sub key {
        my $str = uc(shift);
        $str =~ s|(0*)(\d+)|pack("AN/A*N", "0", $2, length($1))|eg;
        return $str;
    }

    my $heap = Heap::Simple->new(order => "lt",
                                 elements => [Function => \&key]);
    $heap->insert(qw(Athens5.gr Athens40.gr
                     Amsterdam51.nl Amsterdam5.nl amsterdam20.nl));
    print $heap->extract_top, "\n" for 1..$heap->count;
    # This will print:
    Amsterdam5.nl
    amsterdam20.nl
    Amsterdam51.nl
    Athens5.gr
    Athens40.gr

=back

=item X<new_elements>elements => $element_type

This option describes what sort of elements we will store in the heap.
The only reason the module needs to know this is to determine how to access
the key values.

The $element_type is usually an array reference, but if the array has only
one entry, you may use that directly. So you can use:

    elements => "Scalar"

instead of:

    elements => ["Scalar"]

The following element types are currently supported:

=over 2

=item X<Scalar>["Scalar"]

Indicates that the elements are the keys themselves. This is the default if no
elements option is provided. So the constructor in the previous example could
also have been written as:

    my $heap = Heap::Simple->new(order => "lt",
                                 elements => ["Scalar"]);

or in the simplified notation:

    my $heap = Heap::Simple->new(order => "lt", elements => "Scalar");

This element type used to be called C<Key>, and that name is still accepted
for backward compatibility.

=item X<Array>[Array => $index]

Indicates that the elements are array references, with the key at index $index.
So now the element can be not just the key, but also associated data. We can
use this to for example print the values of a hash ordered by key:

    use Heap::Simple;

    my $heap = Heap::Simple->new(order => "lt",
                                 elements => [Array => 0]);
    while (my ($key, $val) = each %hash) {
        $heap->insert([$key, $val]);
    }
    for (1..$heap->count) {
        print $heap->extract_top->[1], "\n";
    }

You can always use something like [$key, @data] to pair up keys and data,
so the "Array" element type is rather generally useful (but see the
L<Object|"Object"> and L<Any|"Any"> element types for another way to pair
keys with data). Since it's so common to have the key in the first position,
you may in fact drop the index in that case, so the constructor in the
previous example could also be written as:

    my $heap = Heap::Simple->new(order => "lt",
                                 elements => ["Array"]);

or using the one element rule:

    my $heap = Heap::Simple->new(order => "lt",
                                 elements => "Array");

In case the elements you want to store are arrays (or array based objects
(or L<fields based objects|fields>) and you are prepared to break the object
encapsulation), this element type is also very nice. If for example the value
on which you want to order is a number at position 4, you could use:

    my $heap = Heap::Simple->new(elements => [Array => 4]);
    print "The key is $object->[4]\n";
    $heap->insert($object);

=item X<Hash>[Hash => $key_name]

Indicates that the elements are hash references, where the key (for the heap
element) is the value $element->{$key_name} .

Redoing the Array example in Hash style gives:

    use Heap::Simple;

    my $heap = Heap::Simple->new(order => "lt",
                                 elements => [Hash => "tag"]);
    while (my ($key, $val) = each %hash) {
        $heap->insert({tag => $key, value => $val});
    }
    for (1..$heap->count) {
        print $heap->extract_top->{value}, "\n";
    }

In case the elements you want to store are hashes (or hash based objects and
you are prepared to break the object encapsulation), this element type is also
very nice. If for example the value on which you want to order is a number
with key "price", you could use:

    my $heap = Heap::Simple->new(elements => [Hash => "price"]);
    print "The key is $object->{price}\n";
    $heap->insert($object);

=item X<Method>[Method => $method_name]

In case you don't want to (or can't) break the object encapsulation, but there
is a method that will return the key for a given object, you can use this.

The method method_name will be called like:

    $key = $element->$method_name();

and should return the key corresponding to $element.

Suppose that the elements are objects whose weight you can access using the
"weight" method. A heap ordered on weight then becomes:

    my $heap = Heap::Simple->new(elements => [Method => "weight"]);
    print "The key is ", $object->weight(), "\n";
    $heap->insert($object);

=item X<Object>[Object => $method_name]

The drawback of the L<"Method" element type|"Method"> is that the method will
be called every time the internals need ordering information, which will be
O(log n) for a single insert or extract on a heap of size n. So it's usually
better to first extract the key before insert, wrap the object with the key
such that key access is cheap and insert that. Since this is so common,
this element type is provided for that.

So this element type will only call $method_name once on the initial insert,
after which internally the key is stored together with the value. This makes
it faster, but it also uses more memory.

It also means that it's now perfectly fine to make changes to the object
that change the key while it is in the heap. This will have absolutely no
influence on the ordering anymore, and methods like L<first_key|"first_key">
will still return what the key value was at insert time.

Repeating the previous example in this style is a trivial variation:

    my $heap = Heap::Simple->new(elements => [Object => "weight"]);
    print "The key is ", $object->weight(), "\n";
    $heap->insert($object);

Since for this element type the key is almost completely decoupled from the
value and only fetched on insert, it often makes sense to not let the heap
calculate the key, but do it yourself before the insert, and then use
L<key_insert|"key_insert">. In fact, if you never use plain L<insert|"insert">
at all, you don't even have to bother passing a method name (though in that
case the fact that the thing you store is an object is pretty irrelevant and
it's probable more natural to use the L<Any element type|"Any">).

=item X<Function>[Function => $code_reference]

For completely general key calculation you can use this element type. The given
code reference will be called on an element like:

    $key = $code_reference->($element);

and should return the key corresponding to $element.

An example:

    sub price {
        my $items = shift;
        my $price = 0;
        $price += $_->price for @$items;
	return $price;
    }

    my $heap = Heap::Simple->new(elements => [Function => \&price]);
    print "All items together will cost ", $item_list->price, "\n";
    $heap->insert($item_list);

=item X<Any>[Any => $code_reference]

The same discussion as under L<Object|"Object"> applies for
L<Function|"Function">: single insert and extract on a size n heap will call
the code reference O(log n) times, which could get slow.

So if you are prepared to use more memory, you can again tell Heap::Simple
to calculate the key already at insert time, and store it together with the
value. This will avoid the need for repeated key calculations.

The "Any" element type will do this for you transparantly.

The heap part of the above example becomes:

    my $heap = Heap::Simple->new(elements => [Any => \&price]);
    print "All items together will cost ", $item_list->price, "\n";
    $heap->insert($item_list);

Since for this element type the key is almost completely decoupled from the
value and only fetched on insert, it often makes sense to not let the heap
calculate the key, but do it yourself before the insert, and then use
L<key_insert|"key_insert">. In fact, if you never use plain L<insert|"insert">
at all, you don't even have to bother passing the code reference. So the last
example could look like:

    my $heap = Heap::Simple->new(elements => "Any");
    my $price = $item_list->price;
    print "All items together will cost $price\n";
    $heap->key_insert($price, $item_list);

Or we can use it to simplify the hash sort on key example a bit:

    use Heap::Simple;

    my $heap = Heap::Simple->new(order => "lt",
                                 elements => "Any");
    # A hash in list context returns a sequence of key/value pairs
    $heap->key_insert(%hash);
    for (1..$heap->count) {
        print $heap->extract_top, "\n";
    }

=back

=item X<new_max_count>max_count => $natural_number

Normally a heap can contain any number of elements. By passing a positive
integer as argument to max_count, you tell the heap that it should never
contain more values than that. If you try to insert something new when the
maximum is reached, the lowest element (with respect to the current order) is
dropped (the thing just being inserted is among the candidates for dropping).

A max count of 0 may or may not be supported depending on the implementor.

You can for example use this to efficiently determine the three highest values
in an array:

    use Heap::Simple;

    my @array = qw(19 3 7 -5 3 18 1);

    my $heap = Heap::Simple->new(max_count => 3);
    $heap->insert(@array);
    print "The three highest values are: ", join(", " => $heap->values), "\n";

    # Will print: The three highest values are: 7, 19, 18

=item X<new_can_die>can_die => $bool

If you use magic values, overload, L<order functions|"CODE"> or
L<key access functions|"new_elements">, then it's possible for external perl
code to run when you do heap operations like L<insert|"insert"> or
L<extract_top|"extract_top">. If these throw an exception, the heap can
be left in an incorrect half changed state.

If you give a true value to can_die, the code for single element operations
will be changed so that they will properly recover by undoing what just got
changed (so a failing operation becomes a no-op). This however will slow down
these operations somewhat, so the default is actually false (most of the time
getting exceptions during the heap operations is impossible anyways).

Operations that insert or extract multiple elements will also get their code
changed so the heap is always left in a consistent state, but the operation
is not atomic since it could already be executed on some of the elements.
You could even lose elements if for example an L<extract_all|"extract_all">
fails halfway through (the already extracted part is gone from the heap but you
never got a chance to store the methods return values).

Multi element operations can be substantially more efficient without this flag
since it may allow the use of better algorithms.

This is a per heap option, so only those heaps that actually set this will
see any slowdown.

All operations that don't change the heap (like L<count|"count"> or
L<top|"top">) are always safe.

Note that all change operations always assume you won't recursively cause
another change to the same heap while they are running. If you do that, all
bets for consistency are are off, even if you set this option.

=item X<dirty>dirty => $bool

Giving this option a true value means that the code may make optimistic
assumptions to gain more speed. These can be things like for example ignoring
overloads, casting all numbers to doubles, ignoring locale for string order,
caching keys etc. The particular optimizations done under dirty will be safe
for most applications though. See the documentation for a particular
implementor (like L<Heap::Simple::XS|Heap::Simple::XS> or
L<Heap::Simple::Perl|Heap::Simple::Perl>) for what the dirty option does for
that package.

The default is no dirty optimizations.

=item X<new_user_data>user_data => $user_data

You can associate one scalar worth of user data with any heap. That scalar can
of course also be a reference to a more complex data structure, allowing you to
effectively associate any amount of data with the heap. This option allows you
to set this scalar value already at object creation. You can use the
L<user_data|"user_data"> method to get/set the associated value.

If this option is not given, the heap starts with "undef" associated to it.

    my $heap = Heap::Simple->new(user_data => "foo");
    print $heap->user_data, "\n";
    # prints foo

=item X<new_infinity>infinity => $infinity

Associates $infinity as the highest possible key with the created heap.
($infinity may or may not be a possible key itself).
Setting it to "undef" means there is no infinity associated with the heap.

The default value depends on the L<order|"new_order"> relation that was
specified.

Usually you can just forget about this option. Only L<top_key|"top_key">
really cares.

=back

Notice that the class into which the resulting heap is blessed will B<not>
be Heap::Simple. It will be an on demand generated class that will have
Heap::Simple as an ancestor.

=item X<insert>$heap->insert(@new_elements)

Inserts each of the @new_elements in the heap. On extraction you get back
exactly the same $element as you inserted, including a possible
L<blessing|perlfunc/"bless">.

In case an exception is raised during insert the heap is only guaranteed to be
in a consistent state if you had set the L<can_die|"new_can_die"> flag to
L<new|"new">. Even then it's possible that some first part of @new_elements has
been inserted into the heap while the rest hasn't (they get inserted in the
order given). You could check how many by calling L<the count method|"count">
before and after the insert. So even with L<can_die|"new_can_die"> only
inserts of single elements are atomic.

Mass insert can be substantially faster if the L<can_die|"new_can_die"> flag
isn't set though.

=item X<key_insert>$heap->key_insert($key1, $element1, $key2, $element2, ...)

Inserts each $element in the heap ordered by the $key given just before it.
Since in this case the key must be stored seperately from the element, this
method only exists for L<"Object"|"Object"> and L<"Any"|"Any"> heaps.

On extraction you get back exactly the same $element as you inserted,
including a possible L<blessing|perlfunc/"bless">.

In case an exception is raised during insert the heap is only guaranteed to be
in a consistent state if you had set the L<can_die|"new_can_die"> flag to
L<new|"new">. Even then it's possible that some first part of the argument list
has been inserted into the heap while the rest hasn't (they get inserted in
the order given). You could check how many by calling
L<the count method|"count"> before and after the insert. So even with
L<can_die|"new_can_die"> only inserts of single key/element pairs are atomic.

Mass insert can be substantially faster if the L<can_die|"new_can_die"> flag
isn't set though.

=item X<extract_top>$element = $heap->extract_top

For all elements in the heap, find the top one (the one that is "lowest" in the
order relation), remove it from the heap and return it.

This method used to be called C<"extract_min"> instead of C<"extract_top">.
The old name is still supported but is deprecated.

Throws an exception if the heap is empty.

=item X<extract_first>$element = $heap->extract_first

Like L<extract_top|"extract_top">, but if the heap is empty it will
return undef (in scalar context).

=item X<top>$element = $heap->top

For all elements in the heap, find the top one (the one that is "lowest" in the
order relation) and return it (without removing it from the heap)..

Throws an exception if the heap is empty.

=item X<first>$element = $heap->first

For all elements in the heap, find the top one (the one that is "lowest" in the
order relation) and return it (without removing it from the heap).For all elements in the heap, find the one with the lowest key and return it.
Returns undef (in scalar context) in case the heap is empty. The contents of
the heap remain unchanged.

Since the data returned from a non-empty heap can often not be undef, you
could use this method to check if a heap is empty, but it's probably more
natural to use L<count|"count"> for that.

Example:

    use Heap::Simple;

    my $heap = Heap::Simple->new;
    $heap->insert(8, 3, 14, -1, 3);
    print $heap->first, "\n";
    # prints -1

=item X<first_key>$top_key = $heap->first_key

Looks for the lowest key in the heap and returns its value. Returns undef
(in scalar context) in case the heap is empty

Example:

    use Heap::Simple;

    my $heap = Heap::Simple->new;
    $heap->insert(8, 3, 14, -1, 3);
    print $heap->first_key, "\n";
    # prints -1

=item X<top_key>$top_key = $heap->top_key

Looks for the lowest key in the heap and returns its value. Returns the highest
possible value (the infinity for the chosen order) in case the heap is empty.
If there is no infinity, it will throw an exception.

Example:

    use Heap::Simple;

    my $heap = Heap::Simple->new;
    $heap->insert(8, 3, 14, -1, 3);
    print $heap->top_key, "\n";
    # prints -1

This method used to be called "min_key" instead of "top_key". The old name is
still supported but is deprecated.

=item X<extract_upto>@elements = $heap->extract_upto($max_key)

Finds all elements in the heap whose key is not above $value and removes them
from the heap (so elements with key equal to $max_key get extracted too).
The list of removed elements is returned ordered by key value (low to high
with repect to the heap order).

Returns an empty list for the empty heap.

Example:

    use Heap::Simple;

    my $heap = Heap::Simple->new;
    $heap->insert(8, 3, 14, -1, 3);
    print join(", ", $heap->extract_upto(3)), "\n";
    # prints -1, 3, 3

This method will lose values in case of an exception even if
L<can_die|"new_can_die"> is true (remember that exceptions of this type are
only possible if you have a self coded key fetch or compare that can die, so
this is normally irrelevant).

=item X<extract_all>@elements = $heap->extract_all

Extracts all elements from $heap and returns them ordered by key value (low to
high with repect to the heap order).

Example:

    use Heap::Simple;

    my $heap = Heap::Simple->new;
    $heap->insert(8, 3, 14, -1, 3);
    print join(", ", $heap->extract_all), "\n";
    # prints -1, 3, 3, 8, 14

This method can lose values in case of an exception even if
L<can_die|"new_can_die"> is true (remember that exceptions of this type are
only possible if you have a key fetch or compare that can die, so this is
normally irrelevant).

If you don't actually care about the order of the elements it's more efficient
to use L<values|"values"> followed by L<clear|"clear">.

It's unspecified what this method returns in scalar context.

=item X<clear>$heap->clear

Removes all elements from the heap. This can be much more efficient than using
L<extract_all|"extract_all"> in a void context.

=item X<count>$count = $heap->count

Returns the number of elements in the heap.
This is an constant time operation, it doesn't really need to count anything.

    use Heap::Simple;

    my $heap = Heap::Simple->new;
    $heap->insert(8, 3, 14, -1, 3);
    print $heap->count, "\n";
    # prints 5

=item X<keys>@keys = $heap->keys

Returns the keys of all elements in the heap in some heap order.
This means that the element at index n is not bigger (in the heap order)
than the element at index 2*n+1 and 2*n+2. So the top key will in fact be
in the first position, but don't expect the whole list to be ordered.

This method may imply a lot of function calls if getting the key from an
element implies a function call (as it does for the L<Method|"Method"> and
L<Function|"Function"> element types, but not for the L<Object|"Object"> and
L<Any|"Any"> element types).

Multiple calls to an unchanged heap will return the keys in the same order,
which is also consistent with the order of L<values|"values">

=item X<values>@values = $heap->values

Returns all elements in the heap in some heap order with respect to the
corresponding keys (see L<keys|"keys">).
Does not remove the values from the heap.

Multiple calls to an unchanged heap will return the values in the same order,
which is also consistent with the order of L<keys|"keys">

=item X<key>$key = $heap->key($value)

Calculates the key corresponding to $value in the same way as the internals
of $heap would. Can fail for L<Object|"Object"> and L<Any|"Any"> element
types if there was no method or function given on heap creation.

Notice that this does not access the elements in the heap in any way.
In particular, it's B<not> looking for $value in the heap hoping to match its
key.

=item X<user_data>$user_data = $heap->user_data

Queries the L<user_data|"new_user_data"> associated with the heap.

=item $old_data = $heap->user_data($new_data)

Associates new L<user_data|"new_user_data"> with the heap. Returns the old
value.

=item X<infinity>$infinity = $heap->infinity

Queries the infinity value associated with the heap. Returns undef if there
is none. The default infinity is implied by the chosen order relation.

=item $old_infinity = $heap->infinity($new_infinity)

Associates a new infinity with the heap. Returns the old value.

=item X<key_index>$key_index = $heap->key_index

Returns the index of the key for L<array reference based heaps|"Array">.
Doesn't exist for the other heap types.

=item X<key_name>$key_name = $heap->key_name

Returns the name of the key key for L<hash reference based heaps|"Hash">.
Doesn't exist for the other heap types.

=item X<key_method>$key_name = $heap->key_method

Returns the name of the method to fetch the key from an object. Only exists
for L<Method|"Method"> and L<Object|"Object"> based heaps.

=item X<key_function>$key_function = $heap->key_function

Returns the code reference of the function to fetch the key from an element.
Only exists for L<"Function"|"Function"> and L<"Any"|"Any"> heaps.

=item X<wrapped>$wrapped => $heap->wrapped

Returns true if key and value are stored seperately (wrapped together in
some internal container), nothing otherwise. This is the sufficient and
necessary condition for L<key_insert|"key_insert"> to work, and will normally
only be true for the L<"Any"|"Any"> and L<Object|"Object"> type heaps.

=item X<max_count>$max_count => $heap->max_count

Returns the maximum size of the heap, or infinity if there is no maximum (the
default unless you used the L<max_count|"new_max_count"> option).

=item X<can_die>$can_die = $heap->can_die

Returns the L<can_die|"new_can_die"> setting for this heap.

=item X<dirty>$dirty = $heap->dirty

Returns the L<dirty|"new_dirty"> setting for this heap.

=item X<order>$order = $heap->order

Returns the L<order|"new_order"> setting for this heap (either C<"E<lt>">, C<"E<gt>">, C<"lt">, C<"gt"> or a code reference).

=item X<elements>@elements = $heap->elements

Returns the L<elements|"new_elements"> setting for this heap. The first entry
in the returned list is a string representing the type in canonical form
(L<"Scalar"|"Scalar"> or L<"Array"|"Array"> etc) followed by any arguments
that type needed (e.g. the key name for a L<Hash|"Hash"> type).

=item X<elements>$elements = $heap->elements

Like the list context version, but only returns the first entry (the canonical
type name).

=item X<absorb>$heap->absorb(@heaps)

Takes all elements from each heap in @heaps and inserts them in $heap,
leaving each heap in @heaps empty. Behaves a bit like:

    for my $work_heap (@heaps) {
        $heap->insert(reverse $work_heap->values);
        $work_heap->clear;
    }

except that it may be more efficient.

If an exception is possible and gets raised during insert, the heaps will be
left in a consistent state with a partial transfer completed on the condition
that L<can_die|"new_can_die"> is set for $heap (the settings for the heaps in
@heaps are irrelevant, their accesses will always be done in a safe way)

=item X<key_absorb>$heap->key_absorb(@heaps)

Takes all elements from each heap in @heaps and key_inserts them in $heap,
leaving each heap in @heaps empty. Behaves a bit like:

    for my $work_heap (@heaps) {
        my @values = $work_heap->values;
        my @keys   = $work_heap->keys;
        $heap->key_insert(pop @keys, pop @values) while @values;
        $work_heap->clear;
    }

except that it's may be more efficient. This is mainly meant for transfer
between wrapped heap types (L<Any|"Any"> and L<Object|"Object">) since it
avoids key recalculation. $heap must of course be a wrapped heap type.

If an exception is possible and gets raised during insert, all heaps will be
left in a consistent state with a partial transfer completed on the condition
that L<can_die|"new_can_die"> is set for $heap (the setting for the heaps in
@heaps are irrelevant, their accesses will always be done in a safe way)

=item X<merge_arrays>my $merged_aref = $heap->merge_arrays($aref1, $aref2, ...)

This convenience function merges a sequence of references to already sorted
arrays into a new sorted array and returns its reference.
So it does something like

    $merge_aref = [sort { $heap->compare_function->($a, $b) } map @$_, @_;
    shift @$merge_aref until @$merge_aref <= $heap->max_count;

except that it's more efficient (e.g. it uses the knowledge that the
argument arrays are already sorted).

It leaves values stored in the $heap completely untouched. $heap is
only used for its attributes: how to find the key, what the compare function is
and the maximum number of elements.

=item X<implementation>$implementation = Heap::Simple->implementation

Returns the package that does the actual work. That will probably be
C<"Heap::Simple::XS"> or C<"Heap::Simple::Perl">.

=back

=head1 SEE ALSO

L<Heap::Simple::Perl>,
L<Heap::Simple::XS>

Some other heap or heap-like classes that exist:

L<Heap>,
L<Heap::Priority>,
L<Array::Heap2>

=head1 AUTHOR

Ton Hospel, E<lt>Heap-Simple@ton.iguana.beE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright 2003 by Ton Hospel

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut
