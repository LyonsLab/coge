package CoGe::Requests::RequestFactory;

sub get {
    my ($self, %message) = @_;

    given ($message->{type}) {
        default       { return; }
    }
}

1;
